package com.hartwig.hmftools.isofox.unmapped;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ConsensusSequence
{
    private final char[] mSequence;
    private int mCurrentLength;

    private final List<UnmappedRead> mCandidateReads;
    private final List<UnmappedRead> mExactMatchReads;
    private final List<UnmappedRead> mInexactMatchReads; // shorter or differing in low-qual bases

    private double mAssignedInexact;

    private final Map<Integer,Map<Character,Integer>> mLowBaseQualCounts; // map of index to base character to count of that base

    private static final int MAX_SEQUENCE = 152;
    private static final int MIN_QUAL = 30;
    private static final char UNSET = '.';

    public ConsensusSequence(final UnmappedRead umRead)
    {
        mSequence = new char[MAX_SEQUENCE];
        mCurrentLength = 0;

        for(int i = 0; i < MAX_SEQUENCE; ++i)
        {
            mSequence[i] = UNSET;
        }

        mLowBaseQualCounts = Maps.newLinkedHashMap();

        mCandidateReads = Lists.newArrayList();
        mExactMatchReads = Lists.newArrayList();
        mInexactMatchReads = Lists.newArrayList();

        initialise(umRead);
    }

    public List<UnmappedRead> exactMatchReads() { return mExactMatchReads; }
    public List<UnmappedRead> inexactMatchReads() { return mInexactMatchReads; }

    public List<UnmappedRead> allMatchedReads()
    {
        final List<UnmappedRead> allReads = Lists.newArrayList(mExactMatchReads);
        allReads.addAll(mInexactMatchReads);
        return allReads;
    }

    public double assignedReadTotal() { return mExactMatchReads.size() + mAssignedInexact; }

    public void addInexactRead(final UnmappedRead umRead, double portion)
    {
        mInexactMatchReads.add(umRead);
        mAssignedInexact += portion;
    }

    public String getSequenceString()
    {
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < mCurrentLength; ++i)
        {
            sb.append(mSequence[i]);
        }

        return sb.toString();
    }

    private void initialise(final UnmappedRead umRead)
    {
        // set details from the first read
        final String bases = umRead.ScBases;
        final byte[] baseQuals = umRead.ScBasesQuals;

        mCurrentLength = bases.length();

        for(int i = 0; i < mCurrentLength; ++i)
        {
            if(baseQuals[i] >= MIN_QUAL)
                setHighQualBase(i, bases.charAt(i));
            else
                addLowQualBase(i, bases.charAt(i), baseQuals[i]);
        }

        mCandidateReads.add(umRead);
    }

    public boolean matches(final UnmappedRead umRead, boolean acceptNewHighQual)
    {
        for(int i = 0; i < min(mCurrentLength, umRead.ScBases.length()); ++i)
        {
            if(umRead.ScBases.charAt(i) == mSequence[i])
                continue;

            if(umRead.ScBasesQuals[i] < MIN_QUAL)
                continue;

            if(isHighQualBase(i))
                return false;

            if(!acceptNewHighQual)
                return false;
        }

        return true;
    }

    public int getBaseDiffCount(final String sequence)
    {
        int diffCount = 0;

        for(int i = 0; i < min(mCurrentLength, sequence.length()); ++i)
        {
            if(sequence.charAt(i) != mSequence[i])
                ++diffCount;
        }

        return diffCount;
    }

    public boolean checkAddReadBases(final UnmappedRead umRead)
    {
        if(mCurrentLength == 0)
        {
            initialise(umRead);
            return true;
        }

        final String bases = umRead.ScBases;
        final byte[] baseQuals = umRead.ScBasesQuals;

        int baseLen = bases.length();

        // initially just check high-qual mismatches
        if(!matches(umRead, true))
            return false;

        /*
        // int highQualMismatches = 0;
        for(int i = 0; i < min(mCurrentLength, baseLen); ++i)
        {
            if(bases.charAt(i) != mSequence[i] && baseQuals[i] >= MIN_QUAL && isHighQualBase(i))
                return false;
        }

        if(highQualMismatches > 0)
            return false;
        */

        // apply the read and cache its details
        mCandidateReads.add(umRead);

        for(int i = 0; i < min(mCurrentLength, baseLen); ++i)
        {
            if(baseQuals[i] >= MIN_QUAL && !isHighQualBase(i))
            {
                // override with any base not previously high-qual
                setHighQualBase(i, bases.charAt(i));
            }
        }

        int excessLength = max(baseLen - mCurrentLength, 0);

        if(excessLength > 0)
        {
            int maxHasHighQualIndex = 0;

            for(int i = mCurrentLength; i < baseLen; ++i)
            {
                if(baseQuals[i] >= MIN_QUAL)
                    maxHasHighQualIndex = i;
            }

            if(maxHasHighQualIndex > 0)
            {
                for(int i = mCurrentLength; i < maxHasHighQualIndex; ++i)
                {
                    if(baseQuals[i] >= MIN_QUAL)
                        setHighQualBase(i, bases.charAt(i));
                    else
                        addLowQualBase(i, bases.charAt(i), baseQuals[i]);
                }

                mCurrentLength = maxHasHighQualIndex;
            }
        }

        return true;
    }

    private boolean isHighQualBase(int index) { return !mLowBaseQualCounts.containsKey(index); }

    private void setHighQualBase(int index, char base)
    {
        mLowBaseQualCounts.remove(index);
        mSequence[index] = base;
    }

    private void setSequenceByQual()
    {
        for(int i = 0; i < mCurrentLength; ++i)
        {
            Map<Character,Integer> baseCounts = mLowBaseQualCounts.get(i);

            if(baseCounts == null)
                continue;

            char maxBase = UNSET;
            int maxQualTotal = 0;

            for(Map.Entry<Character,Integer> entry : baseCounts.entrySet())
            {
                if(maxQualTotal == 0 || entry.getValue() > maxQualTotal)
                {
                    maxQualTotal = entry.getValue();
                    maxBase = entry.getKey();
                }
            }

            mSequence[i] = maxBase;
        }
    }

    private void addLowQualBase(int index, char base, int qual)
    {
        if(mSequence[index] == UNSET)
            mSequence[index] = base;

        Map<Character,Integer> baseCounts = mLowBaseQualCounts.get(index);

        if(baseCounts == null)
        {
            baseCounts = Maps.newHashMap();
            baseCounts.put(base, 1);
            mLowBaseQualCounts.put(index, baseCounts);
            return;
        }

        Integer qualTotal = baseCounts.get(base);
        baseCounts.put(base, qualTotal != null ? qualTotal + qual : qual);
    }

    public double calcHighQualPercent()
    {
        int hqCount = 0;

        for(int i = 0; i < mCurrentLength; ++i)
        {
            if(isHighQualBase(i))
                ++hqCount;
        }

        return hqCount / (double)mCurrentLength;
    }

    public void setExactMatchReads()
    {
        setSequenceByQual();

        for(UnmappedRead umRead : mCandidateReads)
        {
            if(umRead.ScBases.length() < mCurrentLength)
                continue;

            boolean matches = true;

            for(int i = 0; i < mCurrentLength; ++i)
            {
                if(umRead.ScBases.charAt(i) != mSequence[i] && umRead.ScBasesQuals[i] >= MIN_QUAL)
                {
                    matches = false;
                    break;
                }
            }

            if(matches)
            {
                mExactMatchReads.add(umRead);
            }
        }

        mExactMatchReads.forEach(x -> mCandidateReads.remove(x));
    }

    public void reassignReads(final ConsensusSequence other)
    {
        // copy all reads from the other sequence into inexact matches
        other.exactMatchReads().stream()
                .filter(x -> !mExactMatchReads.contains(x) && !mInexactMatchReads.contains(x))
                .forEach(x -> addInexactRead(x, 1));

        other.inexactMatchReads().stream()
                .filter(x -> !mExactMatchReads.contains(x) && !mInexactMatchReads.contains(x))
                .forEach(x -> addInexactRead(x, 1));
    }

    public String toString()
    {
        return String.format("reads(%d) len(%d) seq(%s)",
                mCandidateReads.size(), mCurrentLength, getSequenceString());
    }

}
