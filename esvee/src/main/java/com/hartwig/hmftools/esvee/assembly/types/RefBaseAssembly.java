package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.findUnsetBases;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.read.Read;

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

    private final List<SupportRead> mSupport;

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

        int refBaseLength = assembly.refBaseLength();

        // copy the ref bases from the junction assembly starting at the first ref base (ie the junction base itself)
        int nonJunctionReadExtension = mJunction.isForward() ?
                assembly.minAlignedPosition() - extensionRefPosition : extensionRefPosition - assembly.maxAlignedPosition();

        int copyJuncIndexStart, copyJuncIndexEnd;
        int assemblyIndex;
        int refBaseStart, refBaseEnd;

        if(mJunction.isForward())
        {
            copyJuncIndexEnd = assemblyLength - 1;
            copyJuncIndexStart = copyJuncIndexEnd - refBaseLength + 1;
            refBaseStart = extensionRefPosition;
            refBaseEnd = assembly.minAlignedPosition() - 1;
            assemblyIndex = 0;
            mJunctionSequenceIndex = assemblyLength - 1;
        }
        else
        {
            copyJuncIndexStart = 0;
            copyJuncIndexEnd = refBaseLength - 1;
            refBaseStart = assembly.maxAlignedPosition() + 1;
            refBaseEnd = extensionRefPosition;
            assemblyIndex = assembly.junctionIndex();
            mJunctionSequenceIndex = 0;
        }

        int refGenomeBaseLength = refBaseEnd - refBaseStart + 1;
        byte[] refBases = refGenome != null ? refGenome.getBases(mJunction.Chromosome, refBaseStart, refBaseEnd): new byte[refGenomeBaseLength];
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

        // SV_LOGGER.debug("assembly({}) pre-support bases: {}", assembly, new String(mBases));
    }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public Junction junction() { return mJunction; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }
    public int baseLength() { return mBases.length; }

    public List<SupportRead> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }
    public int readMismatchesCount() { return mReadMismatchesCount; }

    public boolean checkAddRead(final Read read, final SupportType supportType, int permittedMismatches, int requiredOverlap)
    {
        int[] startIndices = getReadAssemblyStartIndices(read);

        if(startIndices == null)
            return false;

        int readStartIndex = startIndices[0];
        int assemblyStartIndex = startIndices[1];

        boolean canAddRead = canAddRead(read, readStartIndex, assemblyStartIndex, permittedMismatches,requiredOverlap);

        if(!canAddRead && readStartIndex < read.getBases().length - 20)
        {
            // run a simple sequence search to find the alignment start where prior indels have offset the read's infer assembly index start
            int readTestEndIndex = readStartIndex + 20;
            int length = readTestEndIndex - readStartIndex + 1;

            if(readStartIndex < 0 || readStartIndex >= read.getBases().length || readStartIndex + length > read.getBases().length)
            {
                SV_LOGGER.error("refAssembly({}) invalid indices({} - {}) vs readBases({}) for ref extension read search",
                        toString(), readStartIndex, readTestEndIndex, read.getBases().length);

                System.exit(1);
            }

            String readBases = new String(read.getBases(), readStartIndex, length);
            assemblyStartIndex = new String(mBases).indexOf(readBases);

            if(assemblyStartIndex >= 0)
            {
                canAddRead = canAddRead(read, readStartIndex, assemblyStartIndex, permittedMismatches, requiredOverlap);
            }
        }

        if(!canAddRead)
        {
            // junction mate reads are added as support even if their ref bases don't match
            if(supportType == SupportType.JUNCTION_MATE)
                mSupport.add(new SupportRead(read, supportType, INVALID_INDEX, 0, permittedMismatches + 1));

            return false;
        }

        addRead(read, supportType, readStartIndex, assemblyStartIndex);

        return true;
    }

    private boolean canAddRead(final Read read, int readStartIndex, int assemblyStartIndex, int permittedMismatches, int requiredOverlap)
    {
        int mismatchCount = 0;
        int overlappedBaseCount = 0;
        int assemblyIndex = assemblyStartIndex;

        for(int i = readStartIndex; i < read.getBases().length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= mBases.length)
                break;

            if(mBases[assemblyIndex] == 0)
                continue;

            ++overlappedBaseCount;

            // any unset base (ie unset qual) can be a mismatch
            byte refBaseQual = mBaseQuals[assemblyIndex] == 0 ? (byte)(LOW_BASE_QUAL_THRESHOLD + 1) : mBaseQuals[assemblyIndex];

            if(!basesMatch(
                    read.getBases()[i], mBases[assemblyIndex], read.getBaseQuality()[i], refBaseQual, LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return overlappedBaseCount >= requiredOverlap;
    }

    private int[] getReadAssemblyStartIndices(final Read read)
    {
        // ensure no indel-adjusted unclipped start is used when aligned to ref bases

        int alignmentStart = read.alignmentStart();
        int leftSoftClipOffset = read.leftClipLength();
        // int unclippedStart = read.unclippedStart();

        if(mJunction.isForward())
        {
            if(!positionWithin(alignmentStart, mExtensionRefPosition, mJunction.Position))
            {
                if(alignmentStart >= mJunction.Position)
                    return null;

                int readIndex = mExtensionRefPosition - alignmentStart + leftSoftClipOffset;
                return new int[] {readIndex, 0};
            }

            return new int[] {leftSoftClipOffset, alignmentStart - mExtensionRefPosition};
        }
        else
        {
            if(!positionWithin(alignmentStart, mJunction.Position, mExtensionRefPosition))
            {
                if(alignmentStart >= mExtensionRefPosition)
                    return null;

                // index off the relative start positions
                int readIndex = mJunction.Position - alignmentStart + leftSoftClipOffset;
                return new int[] {readIndex, mJunctionSequenceIndex};
            }

            return new int[] {leftSoftClipOffset, alignmentStart - mJunction.Position};
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

        mSupport.add(new SupportRead(read, supportType, INVALID_INDEX, 0, mismatchCount));
    }

    public int validRefBaseLength()
    {
        int validRefBaseCount = 0;

        // count bases out from the junction until a gap is encountered
        int index = mJunctionSequenceIndex;

        while(index >= 0 && index < mBases.length)
        {
            if(mBases[index] != 0) //  no check on base qual or actual read base support, despite possible gaps
                ++validRefBaseCount;
            else
                break;

            if(mJunction.isForward())
                --index;
            else
                ++index;
        }

        // and cap at the max observed read's position
        int maxSupportedLength = mJunction.isForward() ?
                mJunction.Position - mMinAlignedPosition + 1 : mMaxAlignedPosition - mJunction.Position + 1;

        return min(maxSupportedLength, validRefBaseCount);
    }

    public void checkValidBases()
    {
        List<int[]> emptyBaseRanges = findUnsetBases(mBases);

        if(!emptyBaseRanges.isEmpty())
        {
            SV_LOGGER.debug("refAssembly({}) has empty ranges: {}", toString(), emptyBaseRanges);
        }
    }

    public String toString()
    {
        return format("junc(%s) aligned(%d - %d) extensionRefPos(%d) length(%d) support(%d) mismatches(%d)",
                mJunction, mMinAlignedPosition, mMaxAlignedPosition, mExtensionRefPosition, baseLength(),
                mSupport.size(), mReadMismatchesCount);
    }
}

