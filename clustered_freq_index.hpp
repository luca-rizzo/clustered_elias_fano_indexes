#pragma once

#include "bitvector_collection.hpp"
#include "compact_elias_fano.hpp"
#include "integer_codes.hpp"
#include "global_parameters.hpp"
#include "configuration.hpp"
#include "clustered_sequence.hpp"
#include "util.hpp"

#include <unordered_set>
#include <fstream>
#include <stdexcept>

namespace ds2i
{
    template <typename FreqsSequence>
    class clustered_freq_index
    {
        typedef clustered_sequence<> cs_t;
        typedef partitioned_sequence<> ps_t;

        public:
        clustered_freq_index()
            : m_num_docs(0)
        {}

        class builder {

        public:
            builder(uint64_t num_docs,
                    global_parameters const& params)
                : m_params(params)
                , m_num_docs(num_docs)
                , m_cur_reference_id(0)
                , m_cur_reference_begin(0)
                , m_cur_reference_size(0)
                
                , m_cur_cluster_postings(0)
                , m_docs_sequences_size(0)
                , m_refs_sequences_size(0)
                
                , m_docs_sequences(params)
                , m_freqs_sequences(params)
                , m_refs_sequences(params)
            {}

            void add_posting_list(ds2i::clustered_binary_collection::sequence const& docs,
                                  ds2i::clustered_binary_collection::sequence const& freqs,
                                  uint64_t occurrences,
                                  sequence_t const& partition)
            {
                size_t size = docs.size();
                m_cur_cluster_postings += size;

                check_size(size);
                if (!m_cur_reference_id) {
                    throw std::runtime_error("Reference must be added first");
                }

                auto const& conf = configuration::get();
                task_region(*conf.executor, [&](task_region_handle& trh)
                {
                    trh.run([&]
                    {
                        succinct::bit_vector_builder docs_bits;
                        write_gamma_nonzero(docs_bits, occurrences);
                        if (occurrences > 1) {
                            docs_bits.append_bits(size, ceil_log2(occurrences + 1));
                        }

                        auto const& freqs_positions =
                            cs_t::write(docs_bits,
                                        m_cur_reference_id - 1,
                                        m_cur_reference_begin,
                                        m_cur_reference_size,
                                        docs, m_num_docs, partition,
                                        m_params);
                        m_docs_sequences.append(docs_bits);

                        sequence_t f;
                        f.reserve(size);
                        for (auto p : freqs_positions) {
                            f.push_back(freqs[p]);
                        }

                        succinct::bit_vector_builder freqs_bits;
                        FreqsSequence::write(freqs_bits,
                                             f.begin(),
                                             occurrences + 1,
                                             size, m_params);
                        m_freqs_sequences.append(freqs_bits);
                    });
                });
            }

            void add_reference(size_t reference_size,
                               reference_iterator_t reference_begin)
            {
                if (m_cur_reference_id) {
                    m_cur_cluster_postings = 0;
                    m_docs_sequences_size = m_docs_sequences.size();
                    m_refs_sequences_size = m_refs_sequences.size();
                }

                check_size(reference_size);
                ++m_cur_reference_id;
                m_cur_reference_begin = reference_begin;
                m_cur_reference_size = reference_size;

                succinct::bit_vector_builder ref_bvb;
                uint64_t log_reference_size = ceil_log2(reference_size);
                write_gamma_nonzero(ref_bvb, log_reference_size);
                if (log_reference_size > 1) {
                    ref_bvb.append_bits(reference_size, log_reference_size + 1);
                }
                cs_t::write_reference(ref_bvb, reference_begin,
                                      m_num_docs, reference_size, m_params);
                m_refs_sequences.append(ref_bvb);
            }

            void build(clustered_freq_index& sq)
            {
                sq.m_num_docs = m_num_docs;
                sq.m_params = m_params;

                m_docs_sequences.build(sq.m_docs_sequences);
                m_freqs_sequences.build(sq.m_freqs_sequences);
                m_refs_sequences.build(sq.m_refs_sequences);
            }

            void print_cluster_stat() {

                std::cout << "\tcluster postings: " << m_cur_cluster_postings << "\n";
                uint64_t cluster_bits =
                    m_docs_sequences.size() - m_docs_sequences_size
                  + m_refs_sequences.size() - m_refs_sequences_size; 
                std::cout << "\tcluster bits x posting: "
                          << double(cluster_bits) / m_cur_cluster_postings
                          << std::endl;
            }

        private:
            global_parameters m_params;
            uint64_t m_num_docs;
            uint64_t m_cur_reference_id;
            reference_iterator_t m_cur_reference_begin;
            size_t m_cur_reference_size;

            size_t m_cur_cluster_postings;
            size_t m_docs_sequences_size;
            size_t m_refs_sequences_size;

            uint32_t m_sequence_id;
            bitvector_collection::builder m_docs_sequences;
            bitvector_collection::builder m_freqs_sequences;
            bitvector_collection::builder m_refs_sequences;

            void check_size(size_t size) {
                if (!size) {
                    throw std::invalid_argument("List must be nonempty");
                }
            }
        };

        uint64_t size() const
        {
            return m_docs_sequences.size();
        }

        uint64_t num_docs() const
        {
            return m_num_docs;
        }
        
        uint64_t clusters() const
        {
            return m_refs_sequences.size();
        }

        class document_enumerator
        {
            public:
            void reset()
            {
                auto val = m_docs_enum.init();
                m_cur_pos = val.first;
                m_cur_docid = val.second;
            }

            void DS2I_FLATTEN_FUNC next()
            {
                auto val = m_docs_enum.next();
                m_cur_pos = val.first;
                m_cur_docid = val.second;
            }

            void DS2I_FLATTEN_FUNC next_geq(uint64_t lower_bound)
            {
                auto val = m_docs_enum.next_geq(lower_bound);
                m_cur_pos = val.first;
                m_cur_docid = val.second;
            }

            uint64_t docid() const
            {
                return m_cur_docid;
            }

            uint64_t DS2I_FLATTEN_FUNC freq()
            {
                return m_freqs_enum.move(m_cur_pos).second;
            }

            uint64_t position() const
            {
                return m_cur_pos;
            }

            uint64_t size() const
            {
                return m_docs_enum.size();
            }

            cs_t::enumerator const& docs_enum() const
            {
                return m_docs_enum;
            }

            typename FreqsSequence::enumerator const& freqs_enum() const
            {
                return m_freqs_enum;
            }

            private:
            friend class clustered_freq_index;

            document_enumerator(cs_t::enumerator docs_enum,
                                typename FreqsSequence::enumerator freqs_enum)
                : m_docs_enum(docs_enum)
                , m_freqs_enum(freqs_enum)
            {
                reset();
            }

            uint64_t m_cur_pos;
            uint64_t m_cur_docid;
            cs_t::enumerator m_docs_enum;
            typename FreqsSequence::enumerator m_freqs_enum;
        };

        document_enumerator operator[](size_t i) const {

            assert(i < size());
            auto docs_it = m_docs_sequences.get(m_params, i);
            uint64_t occurrences = read_gamma_nonzero(docs_it);
            uint64_t n = 1;
            if (occurrences > 1) {
                n = docs_it.take(ceil_log2(occurrences + 1));
            }

            auto freqs_it = m_freqs_sequences.get(m_params, i);
            typename FreqsSequence::enumerator freqs_enum(m_freqs_sequences.bits(),
                                                          freqs_it.position(),
                                                          occurrences + 1, n,
                                                          m_params);
            uint64_t ref_id = read_gamma(docs_it);
            auto ref_it = m_refs_sequences.get(m_params, ref_id);
            cs_t::enumerator docs_enum(m_docs_sequences.bits(),
                                       docs_it.position(),
                                       m_refs_sequences.bits(),
                                       ref_it.position(),
                                       num_docs(), n, m_params);
            return document_enumerator(docs_enum, freqs_enum);
        }

        void warmup(size_t /* i */) const
        {
            // XXX implement this
        }

        global_parameters const& params() const
        {
            return m_params;
        }

        void export_paper_layout(const char* output_filename) const
        {
            if (!output_filename) {
                throw std::invalid_argument("output filename must not be null");
            }

            std::ofstream out(output_filename, std::ios::binary | std::ios::trunc);
            if (!out) {
                throw std::runtime_error("cannot open paper-layout output file");
            }

            auto write_u8 = [&](uint8_t value) {
                out.write(reinterpret_cast<const char*>(&value), sizeof(value));
            };

            auto write_u32 = [&](uint32_t value) {
                out.write(reinterpret_cast<const char*>(&value), sizeof(value));
            };

            auto write_u64 = [&](uint64_t value) {
                out.write(reinterpret_cast<const char*>(&value), sizeof(value));
            };

            auto write_bitvector_words = [&](succinct::bit_vector const& bv) {
                uint64_t bit_size = bv.size();
                uint64_t words = (bit_size + 63) / 64;
                write_u64(bit_size);
                write_u64(words);
                for (uint64_t w = 0; w < words; ++w) {
                    uint64_t word = bv.get_word(w * 64);
                    write_u64(word);
                }
            };

            auto write_stream_section = [&](bitvector_collection const& stream) {
                write_u64(stream.endpoints_count());
                write_u64(stream.endpoints_universe());
                write_bitvector_words(stream.endpoints_bits());
                write_bitvector_words(stream.bits());
            };

            static const char magic[8] = {'C', 'E', 'F', 'P', 'A', 'P', 'R', '1'};
            out.write(magic, sizeof(magic));

            write_u32(2); // format version (v2: endpoints-first stream layout)
            write_u8(1);  // little-endian marker
            write_u8(m_params.ef_log_sampling0);
            write_u8(m_params.ef_log_sampling1);
            write_u8(m_params.rb_log_rank1_sampling);
            write_u8(m_params.rb_log_sampling1);
            write_u8(m_params.log_partition_size);
            write_u8(0);
            write_u8(0);

            write_u64(m_num_docs);
            write_u64(size());
            write_u64(clusters());

            // Three aligned streams in paper order: documents, frequencies, references.
            write_stream_section(m_docs_sequences);
            write_stream_section(m_freqs_sequences);
            write_stream_section(m_refs_sequences);

            out.flush();
            if (!out) {
                throw std::runtime_error("failed while writing paper-layout file");
            }
        }

        void swap(clustered_freq_index& other)
        {
            std::swap(m_params, other.m_params);
            std::swap(m_num_docs, other.m_num_docs);
            m_docs_sequences.swap(other.m_docs_sequences);
            m_freqs_sequences.swap(other.m_freqs_sequences);
            m_refs_sequences.swap(other.m_refs_sequences);
        }

        template <typename Visitor>
        void map(Visitor& visit)
        {
            visit
                (m_params, "m_params")
                (m_num_docs, "m_num_docs")
                (m_docs_sequences, "m_docs_sequences")
                (m_freqs_sequences, "m_freqs_sequences")
                (m_refs_sequences, "m_refs_sequences")
                ;
        }

    private:
        global_parameters m_params;
        uint64_t m_num_docs;
        bitvector_collection m_docs_sequences;
        bitvector_collection m_freqs_sequences;
        bitvector_collection m_refs_sequences;
    };
}
