#ifndef NTHASH__SIMPLE__ITERATOR_H
#define NTHASH__SIMPLE__ITERATOR_H 1

#include <string>
#include <limits>
#include "nthash.hpp"

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */

class ntHashIteratorSimple
{

public:

    /**
     * Default constructor. Creates an iterator pointing to
     * the end of the iterator range.
    */
    ntHashIteratorSimple():
        m_hVec(NULL),
        m_pos(std::numeric_limits<std::size_t>::max())
    {}

    /**
     * Constructor.
     * @param seq address of DNA sequence to be hashed
     * @param k k-mer size
     * @param h number of hashes
    */
    // ntHashIteratorSimple(const std::string& seq, unsigned h, unsigned k): length of seq
    ntHashIteratorSimple(const char *seq, unsigned h, unsigned L, unsigned k, size_t N):
        m_seq(seq), m_seq_length(N), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_L(L), m_pos(0)
    {   
        last_pos = m_seq_length - m_k + 1; 
        m_w = m_L + m_L - k; // 119 length of substring processed by each time
        m_step = m_w - m_L; //119-100=19
        // fprintf(stderr, "m_w: %u\n", m_w); //each time, process a substring of length m_L without non-ATGC.
        if (m_k > m_seq_length) {
            m_pos = std::numeric_limits<std::size_t>::max();
        }
        init();
    }

    ntHashIteratorSimple(const char *seq, size_t pos, unsigned h, unsigned L, unsigned k, size_t N):
        m_seq(seq), m_seq_length(N), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_L(L), m_pos(pos)
    {   
        last_pos = m_seq_length - m_k + 1; 
        m_w = m_L + m_L - k; // 119 length of substring processed by each time
        m_step = m_w - m_L; //119-100=19
        // fprintf(stderr, "m_w: %u\n", m_w); //each time, process a substring of length m_L without non-ATGC.
        if (m_k > m_seq_length) {
            m_pos = std::numeric_limits<std::size_t>::max();
        }
        init();
    }

    inline bool checkATGC(const char *kmerSeq, unsigned& w) {
        // w = 0;
        // while (*kmerSeq != '{' && *kmerSeq != '}' && *kmerSeq != '\0' && w <= m_w) {
        //     ++w;
        // }
        for (w = 0; w < m_w ; ++w) { //at most m_w
            if (*kmerSeq == '{' || *kmerSeq == '}' || *kmerSeq == '\0' ) {
                break;
            }
            ++kmerSeq;
        }
        return w >= m_L; //at least a L-mer
    }

    inline void calcNextWindow(const char *kmerSeq, unsigned& w) {
        for (w = 0; w <= m_step ; ++w) {
            if (*kmerSeq == '{' || *kmerSeq == '}' || *kmerSeq == '\0' ) {
                return;
            }
            ++kmerSeq;
        }
        // return w > 0; //at least a L-mer
    }

    /** Initialize internal state of iterator */
    inline void init() {
        unsigned w = 0;
        // while (m_pos<m_seq_length-m_k+1 && !NTMC64(m_seq+m_pos, m_k, m_h, m_fhVal, m_rhVal, locN, m_hVec, m_hStn))
        // while (m_pos<m_seq_length-m_k+1 && !NTM64(m_seq+m_pos, m_k, m_h, m_hVec, locN)) {
        while (m_pos < last_pos && !checkATGC(m_seq + m_pos, w)) {
            m_pos += w + 1;
        }
        if (m_pos >= last_pos) {
            m_pos = std::numeric_limits<std::size_t>::max();
        } else {
            #ifdef DEBUG
            char *stbuf = (char*)malloc(sizeof(char)*(m_w+1));
            strncpy(stbuf, m_seq + m_pos, w);
            stbuf[w] = '\0';
            fprintf(stderr, "111: %s\n", stbuf);
            // fprintf(stderr, "2222\n");
            // free(stbuf);
            // fprintf(stderr, "22223333\n");
            fprintf(stderr, "m_pos: %u; w: %u;\n", m_pos, w);
            #endif
            // NTM64(m_seq + (m_pos + w - m_L), m_k, m_h, m_hVec);
            m_pos += w - m_L;
            NTM64(m_seq + m_pos, m_k, m_h, m_hVec);
        }
		// for (unsigned i = 0; i < m_h; ++i) {
		// 	m_hStnArray[i] = m_hStn;
		// }
    }

    /** Advance iterator right to the next valid k-mer */
    inline void next()
    {
        // m_pos += m_L - m_k + 1;
        ++ m_pos; //?
        if (m_pos >= last_pos) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        unsigned w = 0;
        calcNextWindow(m_seq + (m_pos + m_L - 1), w);
        // w >= m_L - m_k + 1
        // if (m_k - 1 + w >= m_L) { ///at least m_L
        if (w > 0) { // at least one: m_L - 1 + 1 = m_L; that is at least m_L
            //shift w + 1?
            #ifdef DEBUG
            char *stbuf = (char*)malloc(sizeof(char)*(m_w+1));
            strncpy(stbuf, m_seq + m_pos, m_L - 1 + w);
            stbuf[m_L - 1 + w] = '\0';
            fprintf(stderr, "shift w: %u\n", w);
            fprintf(stderr, "222: %s\n", stbuf);
            // free(stbuf);
            #endif

            NTM64(m_seq + (m_pos - 1), m_seq + (m_pos - 1 + m_k), m_k, m_h, m_hVec, w);
            m_pos += w - 1;
            #ifdef DEBUG
            fprintf(stderr, "m_pos: %u\n", m_pos);
            // exit(0);
            #endif
            // NTM64(m_seq + (m_pos + m_k - 1 + w - m_L), m_k, m_h, m_hVec);
        } else {
            m_pos += m_L - 1;
            init();
        }
    }

    inline size_t pos() const{
    	return m_pos;
    }

    inline uint64_t getHash() {
        // return m_hVec[0];
        const char *str = m_seq + m_pos;
        std::uint64_t hash = m_k;
        for (std::uint32_t j = 0; j < m_k/4; ) {
            std::uint32_t k;
            memcpy(&k, str, 4);
            k += j++;
            hash ^= k;
            hash *= 171717;
            str += 4;
        }
        return hash;
    }
    /** get pointer to hash values for current k-mer */
    inline const uint64_t* operator*() const
    {
        return m_hVec;
    }

    /** test equality with another iterator */
    inline bool operator==(const ntHashIteratorSimple& it) const
    {
        return m_pos == it.m_pos;
    }

    /** test inequality with another iterator */
    inline bool operator!=(const ntHashIteratorSimple& it) const
    {
        return !(*this == it);
    }

    /** pre-increment operator */
    inline ntHashIteratorSimple& operator++()
    {
        next();
        return *this;
    }

    /** iterator pointing to one past last element */
    inline static const ntHashIteratorSimple end()
    {
        return ntHashIteratorSimple();
    }

    /** destructor */
    ~ntHashIteratorSimple() {
        if(m_hVec!=NULL) {
            delete [] m_hVec;
        }
    }

private:

    /** DNA sequence */
    // std::string m_seq;
    const char *m_seq;
    size_t m_seq_length;
    size_t last_pos;

    /** number of hashes */
    unsigned m_h;

    /** k-mer size */
    unsigned m_k;

    /** L L-mer size */
    unsigned m_L;
    unsigned m_w; //window size processed by each time
    unsigned m_step;

    /** hash values */
    uint64_t *m_hVec;
    
    /** position of current k-mer */
    size_t m_pos;
};

#endif
