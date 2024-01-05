package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.read.ReadUtils.getReadPositionAtReferencePosition;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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
    private final BaseMismatches[] mBaseMismatches;

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
        mBaseMismatches = new BaseMismatches[initialAssemblyLength];

        mSupport = Lists.newArrayList();

        addInitialRead(read);
    }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }
    public int[] baseQualTotals() { return mBaseQualTotals; }

    public List<AssemblySupport> support() { return mSupport; }
    public BaseMismatches[] baseMismatchArray() { return mBaseMismatches; }

    public List<BaseMismatch> baseMismatches()
    {
        List<BaseMismatch> mismatches = Lists.newArrayList();
        for(int i = 0; i < mBaseMismatches.length; ++i)
        {
            if(mBaseMismatches[i] != null)
                Arrays.stream(mBaseMismatches[i].Mismatches).filter(x -> x != null).forEach(x -> mismatches.add(x));
        }

        return mismatches;
    }

    private void addInitialRead(final Read read)
    {
        for(int i = 0; i < mBases.length; ++i)
        {
            mBases[i] = 0;
            mBaseQuals[i] = 0;
            mBaseQualTotals[i] = 0;
        }

        tryAddRead(read, true);
    }

    public boolean tryAddRead(final Read read)
    {
        return tryAddRead(read, false);
    }

    private boolean tryAddRead(final Read read, boolean isInitial)
    {
        int junctionReadIndex = read.getReadIndexAtReferencePosition(mInitialJunction.Position, true);

        int readIndexStart;
        int readIndexEnd;
        int assemblyStartIndex;

        if(mInitialJunction.direction() == Direction.FORWARDS)
        {
            readIndexStart = junctionReadIndex;
            readIndexEnd = read.getLength() - 1;
            assemblyStartIndex = 0;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = junctionReadIndex;
            assemblyStartIndex = read.getUnclippedStart() - mMinAlignedPosition;
        }

        int mismatchCount = 0;
        int assemblyIndex = assemblyStartIndex;

        for(int i = readIndexStart; i <= readIndexEnd; ++i, ++assemblyIndex)
        {
            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(isInitial || mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;
                mBaseQualTotals[assemblyIndex] = qual;
            }
            else
            {
                if(mBases[assemblyIndex] == base)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else if(qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else
                {
                    // register a mismatch
                    ++mismatchCount;

                    BaseMismatches baseMismatches = mBaseMismatches[assemblyIndex];

                    if(baseMismatches == null)
                    {
                        baseMismatches = new BaseMismatches(new BaseMismatch(base, read, qual));
                        mBaseMismatches[assemblyIndex] = baseMismatches;
                    }
                    else
                    {
                        baseMismatches.addMismatch(base, read, qual);
                    }
                }
            }
        }

        mSupport.add(new AssemblySupport(read, assemblyStartIndex, readIndexStart, readIndexEnd, mismatchCount));
        return true;
    }

}
