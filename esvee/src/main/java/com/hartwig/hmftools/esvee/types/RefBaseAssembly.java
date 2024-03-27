package com.hartwig.hmftools.esvee.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.types.SupportType.JUNCTION_MATE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.read.Read;

public class RefBaseAssembly
{
    private final Junction mJunction;
    private final int mExtensionRefPosition;
    private int mMinAlignedPosition; // set from supporting reads
    private int mMaxAlignedPosition;
    private int mJunctionSequenceIndex;

    private byte mBases[];
    private byte mBaseQuals[];

    private int mReadMismatchesCount;

    private final List<AssemblySupport> mSupport;

    public RefBaseAssembly(final JunctionAssembly assembly, final int extensionRefPosition, final RefGenomeInterface refGenome)
    {
        mJunction = assembly.junction();
        mExtensionRefPosition = extensionRefPosition;

        // aligned positions will be set from supporting reads only, not assumed
        mMinAlignedPosition = assembly.junction().Position;
        mMaxAlignedPosition = assembly.junction().Position;

        int assemblyLength = abs(extensionRefPosition - mJunction.Position) + 1;

        mBases = new byte[assemblyLength];
        mBaseQuals = new byte[assemblyLength];

        mSupport = Lists.newArrayList();

        mReadMismatchesCount = 0;

        // copy the ref bases from the junction assembly starting at the first ref base (ie the junction base itself)
        // for now ignoring mismatches even if they're dominant could also lower the qual for those, but better to sort this out prior or properly
        int nonJunctionReadExtension = mJunction.isForward() ?
                assembly.minAlignedPosition() - extensionRefPosition : extensionRefPosition - assembly.maxAlignedPosition();

        int copyJuncIndexStart, copyJuncIndexEnd;
        int assemblyIndex;
        int refBaseStart, refBaseEnd;

        if(mJunction.isForward())
        {
            copyJuncIndexStart = nonJunctionReadExtension;
            copyJuncIndexEnd = assemblyLength - 1;
            refBaseStart = extensionRefPosition;
            refBaseEnd = assembly.minAlignedPosition() - 1;
            assemblyIndex = 0;
            mJunctionSequenceIndex = assemblyLength - 1;
        }
        else
        {
            copyJuncIndexStart = 0;
            copyJuncIndexEnd = assembly.maxAlignedPosition() - mJunction.Position + 1;
            refBaseStart = assembly.maxAlignedPosition() + 1;
            refBaseEnd = extensionRefPosition;
            assemblyIndex = assembly.junctionIndex();
            mJunctionSequenceIndex = 0;
        }

        int refBaseLength = refBaseEnd - refBaseStart + 1;
        byte[] refBases = refGenome != null ? refGenome.getBases(mJunction.Chromosome, refBaseStart, refBaseEnd): new byte[refBaseLength];
        int refBaseIndex = 0;

        for(int i = 0; i < mBases.length; ++i)
        {
            if(i >= copyJuncIndexStart && i <= copyJuncIndexEnd && assemblyIndex < assembly.bases().length)
            {
                mBases[i] = assembly.bases()[assemblyIndex];
                mBaseQuals[i] = assembly.baseQuals()[assemblyIndex];
                ++assemblyIndex;
            }
            else
            {
                if(refBaseIndex < refBases.length)
                    mBases[i] = refBases[refBaseIndex++];
                else
                    mBases[i] = 0;

                mBaseQuals[i] = 0;
            }
        }
    }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }
    public int baseLength() { return mBases.length; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }
    public int readMismatchesCount() { return mReadMismatchesCount; }

    public boolean checkAddRead(final Read read, final SupportType supportType, int permittedMismatches, int requiredOverlap)
    {
        int mismatchCount = 0;
        int overlappedBaseCount = 0;

        int[] startIndices = getReadAssemblyStartIndices(read);

        if(startIndices == null)
            return false;

        int readStartIndex = startIndices[0];
        int assemblyStartIndex = startIndices[1];
        int assemblyIndex = assemblyStartIndex;

        for(int i = readStartIndex; i < read.getBases().length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= mBases.length)
                break;

            if(mBases[assemblyIndex] == 0)
                continue;

            ++overlappedBaseCount;

            if(!basesMatch(
                    read.getBases()[i], mBases[assemblyIndex], read.getBaseQuality()[i], mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches && supportType != JUNCTION_MATE)
                    return false;
            }
        }

        if(overlappedBaseCount < requiredOverlap)
            return false;

        addRead(read, supportType, readStartIndex, assemblyStartIndex);

        return true;
    }

    public int validRefBaseLength()
    {
        int validRefBaseCount = 0;

        // count bases out from the junction until a gap is encountered
        int index = mJunctionSequenceIndex;

        while(index >= 0 && index < mBases.length)
        {
            if(mBases[index] != 0)
                ++validRefBaseCount;
            else
                break;

            if(mJunction.isForward())
                --index;
            else
                ++index;
        }

        return validRefBaseCount;
    }

    private int[] getReadAssemblyStartIndices(final Read read)
    {
        if(mJunction.isForward())
        {
            if(!positionWithin(read.unclippedStart(), mExtensionRefPosition, mJunction.Position))
            {
                if(read.unclippedStart() >= mJunction.Position)
                    return null;

                int readIndex = mExtensionRefPosition - read.unclippedStart();
                return new int[] {readIndex, 0};
            }

            return new int[] {0, read.unclippedStart() - mExtensionRefPosition};
        }
        else
        {
            if(!positionWithin(read.unclippedStart(), mJunction.Position, mExtensionRefPosition))
            {
                if(read.unclippedStart() >= mExtensionRefPosition)
                    return null;

                // index off the relative start positions
                int readIndex = mJunction.Position - read.unclippedStart();
                return new int[] {readIndex, mJunctionSequenceIndex};
            }

            return new int[] {0, read.unclippedStart() - mJunction.Position};
        }
    }

    private void addRead(final Read read, final SupportType supportType, int readStartIndex, int assemblyStartIndex)
    {
        int mismatchCount = 0;

        int[] readIndexRange = {readStartIndex, read.getBases().length - 1};
        int assemblyIndex = assemblyStartIndex;

        for(int i = readStartIndex; i < read.getBases().length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= mBases.length)
            {
                readIndexRange[1] = i;
                break;
            }

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

                    ++mReadMismatchesCount;
                }
            }
        }

        if(mJunction.isForward())
        {
            mMinAlignedPosition = min(mMinAlignedPosition, read.alignmentStart());
        }
        else
        {
            mMaxAlignedPosition = max(mMaxAlignedPosition, read.alignmentEnd());
        }

        mSupport.add(new AssemblySupport(read, supportType, assemblyIndex, 0, 0, mismatchCount));
    }

    public String toString()
    {
        return format("junc(%s) aligned(%d - %d) extensionRefPos(%d) length(%d) support(%d) mismatches(%d)",
                mJunction, mMinAlignedPosition, mMaxAlignedPosition, mExtensionRefPosition, baseLength(),
                mSupport.size(), mReadMismatchesCount);
    }
}

