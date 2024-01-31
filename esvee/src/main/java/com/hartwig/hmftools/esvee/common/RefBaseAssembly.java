package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class RefBaseAssembly
{
    private final Junction mJunction;
    private final int mExtensionRefPosition;
    private final int nNonJunctionReadExtension;
    private int mJunctionSequenceIndex;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<AssemblySupport> mSupport;

    private final SequenceMismatches mSequenceMismatches;

    public RefBaseAssembly(final JunctionAssembly assembly, final int extensionRefPosition)
    {
        mJunction = assembly.junction();
        mExtensionRefPosition = extensionRefPosition;

        int assemblyLength = abs(extensionRefPosition - mJunction.Position) + 1;

        mBases = new byte[assemblyLength];
        mBaseQuals = new byte[assemblyLength];

        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();

        // copy the ref bases from the junction assembly, for now ignoring mismatches even if they're dominant
        // could also lower the qual for those, but better to sort this out prior or properly

        nNonJunctionReadExtension = mJunction.isForward() ?
                assembly.minAlignedPosition() - extensionRefPosition : extensionRefPosition - assembly.maxAlignedPosition();

        int copyIndexStart;
        int copyIndexEnd;
        int assemblyIndex;

        if(mJunction.isForward())
        {
            copyIndexStart = nNonJunctionReadExtension;
            copyIndexEnd = assemblyLength - 1;
            assemblyIndex = 0;
            mJunctionSequenceIndex = assemblyLength - 1;
        }
        else
        {
            copyIndexStart = 0;
            copyIndexEnd = assembly.maxAlignedPosition() - mJunction.Position + 1;
            assemblyIndex = assembly.junctionIndex();
            mJunctionSequenceIndex = 0;
        }

        for(int i = 0; i < mBases.length; ++i)
        {
            if(i >= copyIndexStart && i <= copyIndexEnd && assemblyIndex < assembly.bases().length)
            {
                mBases[i] = assembly.bases()[assemblyIndex];
                mBaseQuals[i] = assembly.baseQuals()[assemblyIndex];
                ++assemblyIndex;
            }
            else
            {
                mBases[i] = 0;
                mBaseQuals[i] = 0;
            }
        }
    }

    public int extensionRefPosition() { return mExtensionRefPosition; }
    public int baseLength() { return mBases.length; }

    public int nonJunctionReadExtension() { return nNonJunctionReadExtension; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public SequenceMismatches mismatches() { return mSequenceMismatches; }

    private int getReadAssemblyStartIndex(final Read read)
    {
        if(mJunction.isForward())
        {
            if(!positionsWithin(read.unclippedStart(), read.unclippedEnd(), mExtensionRefPosition, mJunction.Position))
                return INVALID_INDEX;

            return read.unclippedStart() - mExtensionRefPosition;
        }
        else
        {
            if(!positionsWithin(read.unclippedStart(), read.unclippedEnd(), mJunction.Position, mExtensionRefPosition))
                return INVALID_INDEX;

            return read.unclippedStart() - mJunction.Position;
        }
    }

    private static final int INVALID_INDEX = -1;

    public boolean checkAddRead(final Read read, final SupportType supportType, int permittedMismatches)
    {
        int mismatchCount = 0;

        int assemblyStartIndex = getReadAssemblyStartIndex(read);

        if(assemblyStartIndex == INVALID_INDEX)
            return false;

        int assemblyIndex = assemblyStartIndex;

        for(int i = 0; i < read.getBases().length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= mBases.length)
                break;

            if(!basesMatch(
                    read.getBases()[i], mBases[assemblyIndex], read.getBaseQuality()[i], mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        addRead(read, supportType, assemblyStartIndex);

        return true;
    }

    private void addRead(final Read read, final SupportType supportType, int assemblyStartIndex)
    {
        int mismatchCount = 0;

        int[] readIndexRange = {0, read.getBases().length - 1};
        int assemblyIndex = assemblyStartIndex;

        for(int i = 0; i < read.getBases().length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= mBases.length)
                break;

            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                }
                else
                {
                    ++mismatchCount;
                    mSequenceMismatches.add(assemblyIndex, base, read, qual);
                }
            }
        }

        mSupport.add(new AssemblySupport(read, supportType, assemblyIndex, 0, readIndexRange, mismatchCount));
    }
}
