package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;

import java.util.Arrays;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class AssemblySequence
{
    private final Junction mInitialJunction;
    private final Read mInitialRead;

    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];
    private int mBaseQualTotals[];

    private final List<AssemblySupport> mSupport;

    private final SequenceMismatches mSequenceMismatches;

    public AssemblySequence(final Junction initialJunction, final Read read, final int minAlignedPosition, final int maxAlignedPosition)
    {
        mInitialJunction = initialJunction;
        mInitialRead = read;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        int initialAssemblyLength = maxAlignedPosition - minAlignedPosition + 1;
        mBases = new byte[initialAssemblyLength];

        mBaseQuals = new byte[initialAssemblyLength];
        mBaseQualTotals = new int[initialAssemblyLength];
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();

        addInitialRead(read);
    }

    public Read initialRead() { return mInitialRead; }
    public Junction initialJunction() { return mInitialJunction; }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int length() { return mMaxAlignedPosition - mMinAlignedPosition + 1; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }
    public int[] baseQualTotals() { return mBaseQualTotals; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }
    public SequenceMismatches mismatches() { return mSequenceMismatches; }

    private void addInitialRead(final Read read)
    {
        for(int i = 0; i < mBases.length; ++i)
        {
            mBases[i] = 0;
            mBaseQuals[i] = 0;
            mBaseQualTotals[i] = 0;
        }

        addRead(read, false);
    }

    public boolean checkReadMatches(final Read read, int permittedMismatches)
    {
        ReadSequenceCoords readCoords = new ReadSequenceCoords(read, this);

        int mismatchCount = 0;

        int assemblyIndex = readCoords.AssemblyStartIndex;

        for(int i = readCoords.ReadIndexStart; i <= readCoords.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            if(!basesMatch(
                    read.getBases()[i], read.getBaseQuality()[i], mBases[assemblyIndex], mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return true;
   }

    private class ReadSequenceCoords
    {
        public final int JunctionReadIndex;
        public final int AssemblyStartIndex;
        public final int ReadIndexStart;
        public final int ReadIndexEnd;

        public ReadSequenceCoords(final Read read, final AssemblySequence sequence)
        {
            JunctionReadIndex = read.getReadIndexAtReferencePosition(sequence.initialJunction().Position, true);

            if(mInitialJunction.isForward())
            {
                ReadIndexStart = JunctionReadIndex;
                ReadIndexEnd = read.getLength() - 1;
                AssemblyStartIndex = 0;
            }
            else
            {
                ReadIndexStart = 0;
                ReadIndexEnd = JunctionReadIndex;
                AssemblyStartIndex = read.getUnclippedStart() - mMinAlignedPosition;
            }
        }

        public String toString()
        {
            return format("asmIndex(%d) read(junc=%d start=%d end=%d)",
                AssemblyStartIndex, JunctionReadIndex, ReadIndexStart, ReadIndexEnd);
        }
    }

    public void addRead(final Read read, boolean registerMismatches)
    {
        ReadSequenceCoords readCoords = new ReadSequenceCoords(read, this);

        int mismatchCount = 0;
        int assemblyIndex = readCoords.AssemblyStartIndex;

        for(int i = readCoords.ReadIndexStart; i <= readCoords.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;
                mBaseQualTotals[assemblyIndex] = qual;
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else if(registerMismatches)
                {
                    // register a mismatch
                    ++mismatchCount;

                    if(registerMismatches)
                        mSequenceMismatches.add(assemblyIndex, base, read, qual);
                }
            }
        }

        mSupport.add(new AssemblySupport(
                read, readCoords.AssemblyStartIndex, readCoords.ReadIndexStart, readCoords.ReadIndexEnd, mismatchCount));
    }

    public boolean matches(final AssemblySequence other, int minOverlapBases, int maxMismatchCount)
    {
        int overlapMinIndex = max(mMinAlignedPosition, other.minAlignedPosition());
        int overlapMaxIndex = min(mMaxAlignedPosition, other.maxAlignedPosition());

        int overlapLength = overlapMaxIndex - overlapMinIndex + 1;
        if(overlapLength < minOverlapBases)
            return false;

        int offset = overlapMinIndex - minAlignedPosition();
        int otherOffset = overlapMinIndex - other.minAlignedPosition();

        int mismatchCount = 0;

        for(int i = 0; i < overlapLength; ++i)
        {
            if(mBases[i + offset] != other.bases()[i + otherOffset])
            {
                ++mismatchCount;

                if(mismatchCount > maxMismatchCount)
                    return false;
            }
        }

        return true;
    }


}
