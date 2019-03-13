/* from https://github.com/bcgsc/ntHash
 * some changes are made
 */
/*
    MIT License

    Copyright (c) 2018 Hamid Mohamadi

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

/*
 *
 * ntHashIterator.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */
#ifndef NTHASH__ITERATOR_H
#define NTHASH__ITERATOR_H 1

#include <string>
#include <limits>
#include <omp.h>
#include <atomic>
#include "nthash.hpp"

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */

class ntHashIterator
{

public:

    /**
     * Default constructor. Creates an iterator pointing to
     * the end of the iterator range.
    */
    ntHashIterator():
        m_hVec(NULL),
        m_pos(std::numeric_limits<std::size_t>::max())
    {}

    /**
     * Constructor.
     * @param seq address of DNA sequence to be hashed
     * @param k k-mer size
     * @param h number of hashes
    */
    // ntHashIterator(const std::string& seq, unsigned h, unsigned k): length of seq
    /*ntHashIterator(const char *seq, unsigned h, unsigned k, size_t N):
        m_seq(seq), m_seq_length(N), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_pos(0)
    {   
        if (m_k > m_seq_length) {
            m_pos = std::numeric_limits<std::size_t>::max();
        }
        init();
    }*/ 

    ntHashIterator(const char *seq, size_t pos, unsigned h, unsigned k, size_t N):
        m_seq(seq), m_seq_length(N), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_pos(pos)
    {   
        if (m_k > m_seq_length) {
            m_pos = std::numeric_limits<std::size_t>::max();
        }
        init();
    }

    /** Initialize internal state of iterator */
    inline void init()
    {
        if (m_k > m_seq_length) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        unsigned locN=0;
        // while (m_pos<m_seq.length()-m_k+1 && !NTM64(m_seq.data()+m_pos, m_k, m_h, locN, m_hVec, m_hStn))
        while (m_pos<m_seq_length-m_k+1 && !NTM64(m_seq+m_pos, m_k, m_h, m_hVec, locN))
            m_pos+=locN+1;
        // NTM64(m_seq+m_pos, m_k, m_h, m_hVec);
        /*char *stbuf = (char*)malloc(sizeof(char)*(m_k+1));
        strncpy(stbuf, m_seq + m_pos, m_k);
        stbuf[m_k] = '\0';
        fprintf(stderr, "%s\n", stbuf);*/

        if (m_pos >= m_seq_length-m_k+1)
            m_pos = std::numeric_limits<std::size_t>::max();
    }

    /** Advance iterator right to the next valid k-mer */
    inline void next()
    {
        ++m_pos;
        if (m_pos >= m_seq_length-m_k+1) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        // if(seedTab[(unsigned char)(m_seq[m_pos+m_k-1])]==seedN) { //changes
        if (m_seq[m_pos+m_k-1] == '{' || m_seq[m_pos+m_k-1] == '}' || m_seq[m_pos+m_k-1] == '\0') {
            m_pos += m_k;
            init();
        }
        else { // changes
            /*char *stbuf = (char*)malloc(sizeof(char)*(m_k+1));
            strncpy(stbuf, m_seq + m_pos, m_k);
            stbuf[m_k] = '\0';
            fprintf(stderr, "%s\n", stbuf);*/

            NTM64(m_seq[m_pos-1], m_seq[m_pos-1+m_k], m_k, m_h, m_hVec);
        }
        
    }

    inline size_t pos() const{
    	return m_pos;
    }

    // added by Yuansheng Liu
    inline uint64_t getHash() { 
        // return m_hVec[0];
        // 
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
    inline bool operator==(const ntHashIterator& it) const
    {
        return m_pos == it.m_pos;
    }

    /** test inequality with another iterator */
    inline bool operator!=(const ntHashIterator& it) const
    {
        return !(*this == it);
    }

    /** pre-increment operator */
    inline ntHashIterator& operator++()
    {
        next();
        return *this;
    }

    /** iterator pointing to one past last element */
    inline static const ntHashIterator end()
    {
        return ntHashIterator();
    }

    /** destructor */
    ~ntHashIterator() {
        if(m_hVec!=NULL) {
            delete [] m_hVec;
        }
    }

private:

    /** DNA sequence */
    // std::string m_seq;
    const char *m_seq;
    size_t m_seq_length;

    /** number of hashes */
    unsigned m_h;

    /** k-mer size */
    unsigned m_k;

    /** hash values */
    uint64_t *m_hVec;
    
    /** position of current k-mer */
    size_t m_pos;
};

#endif
