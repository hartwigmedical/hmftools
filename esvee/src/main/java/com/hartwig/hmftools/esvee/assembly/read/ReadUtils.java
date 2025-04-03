package com.hartwig.hmftools.esvee.assembly.read;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.maxConcordantFragmentLength;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.common.CommonUtils;

import htsjdk.samtools.CigarElement;

public final class ReadUtils
{
    public static boolean isDiscordantFragment(final Read read)
    {
        return CommonUtils.isDiscordantFragment(
                read.bamRecord(), maxConcordantFragmentLength(MAX_OBSERVED_CONCORDANT_FRAG_LENGTH), read.supplementaryData());
    }

    public static boolean isValidSupportCoordsVsJunction(final Read read, boolean isForwardJunction, int junctionPosition)
    {
        // cannot cross the junction since will already have considered all junction candidate reads
        // and must read in the direction of the junction, and cannot be soft-clipped prior to the junction
        if(isForwardJunction)
        {
            if(read.negativeStrand() || read.isRightClipped())
                return false;

            if(read.alignmentEnd() > junctionPosition)
                return false;
        }
        else
        {
            if(read.positiveStrand() || read.isLeftClipped())
                return false;

            if(read.alignmentStart() < junctionPosition)
                return false;
        }

        return true;
    }

    public static final int INVALID_INDEX = -1;

    public static int getReadIndexAtReferencePosition(final Read read, final int refPosition, boolean allowExtrapolation)
    {
        return getReadIndexAtReferencePosition(read, refPosition, allowExtrapolation, null);
    }

    public static int getReadIndexAtReferencePosition(final Read read, final Junction junction, boolean allowExtrapolation)
    {
        return getReadIndexAtReferencePosition(read, junction.Position, allowExtrapolation, junction.Orient);
    }

    public static int getReadIndexAtReferencePosition(
            final Read read, final int refPosition, boolean allowExtrapolation, @Nullable final Orientation requiredOrientation)
    {
        // finds the read index given a reference position, and extrapolates outwards from alignments as required

        int alignmentStart = read.alignmentStart();
        int alignmentEnd = read.alignmentEnd();

        // for indel reads, use the implied alignment and unclipped positions for the applicable direction
        if((requiredOrientation == null || requiredOrientation.isReverse()) && allowExtrapolation && read.indelImpliedAlignmentStart() > 0)
            alignmentStart = read.indelImpliedAlignmentStart();

        if((requiredOrientation == null || requiredOrientation.isForward()) && allowExtrapolation && read.indelImpliedAlignmentEnd() > 0)
            alignmentEnd = read.indelImpliedAlignmentEnd();

        if(refPosition <= alignmentStart)
        {
            if(!allowExtrapolation && refPosition < alignmentStart)
                return INVALID_INDEX;

            int baseDiff = alignmentStart - refPosition;
            int unclippedStart = read.indelImpliedAlignmentStart() > 0 ? read.indelImpliedUnclippedStart() : read.unclippedStart();
            int softClipBases = alignmentStart - unclippedStart;
            return baseDiff <= softClipBases ? softClipBases - baseDiff : INVALID_INDEX;
        }
        else if(refPosition >= alignmentEnd)
        {
            if(!allowExtrapolation && refPosition > alignmentEnd)
                return INVALID_INDEX;

            int baseDiff = refPosition - alignmentEnd;
            int unclippedEnd = read.indelImpliedAlignmentEnd() > 0 ? read.indelImpliedUnclippedEnd() : read.unclippedEnd();
            int softClipBases = unclippedEnd - alignmentEnd;
            return baseDiff <= softClipBases ? read.basesLength() - (softClipBases - baseDiff) - 1 : INVALID_INDEX;
        }

        // cannot use standard method since CIGAR and coords may have been adjusted
        int readIndex = 0;
        int currentPos = read.alignmentStart();
        for(CigarElement element : read.cigarElements())
        {
            if(!element.getOperator().consumesReferenceBases())
            {
                if(element.getOperator().consumesReadBases())
                    readIndex += element.getLength();

                continue;
            }

            if(currentPos == refPosition)
                break;

            if(!element.getOperator().consumesReadBases())
            {
                // for a D or N where the position is inside it, return the read index for the start of the element
                if(refPosition >= currentPos && refPosition < currentPos + element.getLength())
                    return readIndex - 1;

                currentPos += element.getLength();
            }
            else
            {
                // pos = 100, element = 10M, covering pos 100-109, read index 4 (say after 4S), ref pos at last base of element = 109
                if(refPosition >= currentPos && refPosition < currentPos + element.getLength())
                    return readIndex + refPosition - currentPos;

                currentPos += element.getLength();
                readIndex += element.getLength();
            }
        }

        return readIndex;
    }

    public static int avgBaseQuality(final Read read) { return avgBaseQuality(read.getBaseQuality(), 0, read.basesLength() - 1); }

    public static int avgBaseQuality(final byte[] baseQualities, final int indexStart, final int indexEnd)
    {
        if(indexStart > indexEnd || indexStart < 0 || indexEnd >= baseQualities.length)
            return -1;

        int qualitySum = 0;
        for(int i = indexStart; i <= indexEnd; i++)
        {
            qualitySum += baseQualities[i];
        }

        return qualitySum / (indexEnd - indexStart + 1);
    }
}
