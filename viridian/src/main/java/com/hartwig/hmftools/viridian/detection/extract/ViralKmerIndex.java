package com.hartwig.hmftools.viridian.detection.extract;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_KMER_BLOOM_BITS;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_KMER_BLOOM_HASHES;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_KMER_LENGTH;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// k-mer based Bloom filter to screen candidate viral reads.
// If a read has a k-mer which is in a virus genome, then it's probably viral, otherwise it's probably not viral.
// This is purely a performance optimisation. Without it, a typical sample will have gigabytes of false positive candidate viral reads.
public class ViralKmerIndex
{
    private static final int INVALID_BASE = -1;
    private static final long NO_KMER = -1L;
    private static final int[] BASE_CODES = buildBaseCodes();

    private final int mKmerLength;
    private final long mKmerMask;
    private final int mTopBaseShift;
    private final long[] mBloomWords;
    private final long mBloomBitMask;
    private long mKmerCount;

    private ViralKmerIndex(int kmerLength, int bloomBits)
    {
        mKmerLength = kmerLength;
        mKmerMask = (1L << (2 * kmerLength)) - 1;
        mTopBaseShift = 2 * (kmerLength - 1);
        mBloomWords = new long[bloomBits >>> 6];
        mBloomBitMask = bloomBits - 1L;
    }

    public long kmerCount() { return mKmerCount; }

    public boolean hasViralKmer(byte[] bases)
    {
        KmerIterator iterator = new KmerIterator(bases);
        long kmer;
        while((kmer = iterator.nextKmer()) != NO_KMER)
        {
            if(bloomContains(kmer))
            {
                return true;
            }
        }
        return false;
    }

    public static ViralKmerIndex build(String virusFastaFile)
    {
        ViralKmerIndex index = new ViralKmerIndex(VIRAL_KMER_LENGTH, VIRAL_KMER_BLOOM_BITS);
        try(IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(virusFastaFile)))
        {
            for(SAMSequenceRecord sequence : fasta.getSequenceDictionary().getSequences())
            {
                index.addSequence(fasta.getSequence(sequence.getSequenceName()).getBases());
            }
            return index;
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to read virus reference FASTA for k-mer index", e);
        }
    }

    static ViralKmerIndex fromSequences(List<byte[]> sequences)
    {
        ViralKmerIndex index = new ViralKmerIndex(VIRAL_KMER_LENGTH, VIRAL_KMER_BLOOM_BITS);
        sequences.forEach(index::addSequence);
        return index;
    }

    private void addSequence(byte[] bases)
    {
        KmerIterator iterator = new KmerIterator(bases);
        long kmer;
        while((kmer = iterator.nextKmer()) != NO_KMER)
        {
            bloomAdd(kmer);
        }
    }

    private final class KmerIterator
    {
        private final byte[] mBases;
        private int mPosition;
        private long mForward;
        private long mReverse;
        private int mValidBasesInWindow;

        private KmerIterator(byte[] bases)
        {
            mBases = bases;
        }

        private long nextKmer()
        {
            while(mPosition < mBases.length)
            {
                int code = BASE_CODES[mBases[mPosition++] & 0xFF];
                if(code == INVALID_BASE)
                {
                    mForward = 0;
                    mReverse = 0;
                    mValidBasesInWindow = 0;
                    continue;
                }
                mForward = ((mForward << 2) | code) & mKmerMask;
                mReverse = (mReverse >>> 2) | ((long) (3 - code) << mTopBaseShift);
                if(++mValidBasesInWindow >= mKmerLength)
                {
                    return Math.min(mForward, mReverse);
                }
            }
            return NO_KMER;
        }
    }

    private boolean bloomContains(long kmer)
    {
        long hash = kmer;
        for(int i = 0; i < VIRAL_KMER_BLOOM_HASHES; ++i)
        {
            hash = mix(hash);
            if(!bitIsSet(hash & mBloomBitMask))
            {
                return false;
            }
        }
        return true;
    }

    private void bloomAdd(long kmer)
    {
        long hash = kmer;
        for(int i = 0; i < VIRAL_KMER_BLOOM_HASHES; ++i)
        {
            hash = mix(hash);
            setBit(hash & mBloomBitMask);
        }
        ++mKmerCount;
    }

    private boolean bitIsSet(long bit)
    {
        return (mBloomWords[(int) (bit >>> 6)] & (1L << (bit & 63))) != 0;
    }

    private void setBit(long bit)
    {
        mBloomWords[(int) (bit >>> 6)] |= 1L << (bit & 63);
    }

    private static long mix(long value)
    {
        value ^= value >>> 33;
        value *= 0xff51afd7ed558ccdL;
        value ^= value >>> 33;
        value *= 0xc4ceb9fe1a85ec53L;
        value ^= value >>> 33;
        return value;
    }

    private static int[] buildBaseCodes()
    {
        int[] codes = new int[256];
        Arrays.fill(codes, INVALID_BASE);
        codes['A'] = codes['a'] = 0;
        codes['C'] = codes['c'] = 1;
        codes['G'] = codes['g'] = 2;
        codes['T'] = codes['t'] = 3;
        return codes;
    }
}
