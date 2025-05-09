package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNC_POS_DIFF;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.maxConcordantFragmentLength;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
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

    public static int getReadIndexAtJunction(final Read read, final Junction junction, boolean allowExtrapolation)
    {
        if(junction.indelBased() && read.indelCoords() != null)
        {
            return readIndexFromPosition(read, junction.Position);
        }

        return getReadIndexAtReferencePosition(read, junction.Position, allowExtrapolation, junction.Orient);
    }

    public static int getReadIndexAtReferencePosition(
            final Read read, final int refPosition, boolean allowExtrapolation, @Nullable final Orientation requiredOrientation)
    {
        // finds the read index given a reference position, and extrapolates outwards from alignments as required

        int alignmentStart = read.alignmentStart();
        int alignmentEnd = read.alignmentEnd();

        // for indel reads, use the implied alignment and unclipped positions for the applicable direction
        if((requiredOrientation == null || requiredOrientation.isReverse()) && read.hasIndelImpliedUnclippedStart())
        {
            if(read.indelImpliedUnclippedStart() > refPosition)
                return INVALID_INDEX;

            return refPosition - read.indelImpliedUnclippedStart();
        }

        if((requiredOrientation == null || requiredOrientation.isForward()) && read.hasIndelImpliedUnclippedEnd())
        {
            if(read.indelImpliedUnclippedEnd() < refPosition)
                return INVALID_INDEX;

            return read.indelImpliedUnclippedEnd() - refPosition - 1;
        }

        // use the unclipped position to determine how far back from the start or end of the read the junction index is
        if(refPosition <= alignmentStart)
        {
            if(!allowExtrapolation && refPosition < alignmentStart)
                return INVALID_INDEX;

            int baseDiff = alignmentStart - refPosition;
            int softClipBases = alignmentStart - read.unclippedStart();
            return baseDiff <= softClipBases ? softClipBases - baseDiff : INVALID_INDEX;
        }
        else if(refPosition >= alignmentEnd)
        {
            if(!allowExtrapolation && refPosition > alignmentEnd)
                return INVALID_INDEX;

            int baseDiff = refPosition - alignmentEnd;
            int softClipBases = read.unclippedEnd() - alignmentEnd;
            return baseDiff <= softClipBases ? read.basesLength() - (softClipBases - baseDiff) - 1 : INVALID_INDEX;
        }

        // cannot use standard method since CIGAR and coords may have been adjusted
        return readIndexFromPosition(read, refPosition);
    }

    private static int readIndexFromPosition(final Read read, final int refPosition)
    {
        // could replace this with a standard method if as efficient
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

    public static boolean readSoftClipsAndCrossesJunction(final Read read, final Junction junction, final RefGenomeInterface refGenome)
    {
        if(junction.isForward())
        {
            // first check indel-inferred soft-clips
            if(read.hasIndelImpliedUnclippedEnd() && read.maxUnclippedEnd() > junction.Position)
                return true;

            // soft-clip must be close enough to the junction
            if(read.isRightClipped())
                return abs(read.alignmentEnd() - junction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF && read.unclippedEnd() > junction.Position;
        }
        else
        {
            if(read.hasIndelImpliedUnclippedStart() && read.minUnclippedStart() < junction.Position)
                return true;

            if(read.isLeftClipped())
                return abs(read.alignmentStart() - junction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF && read.unclippedStart() < junction.Position;
        }

        // check for support from ref bases with mismatches
        return refGenome != null ? readMismatchesPastJunction(read, junction, refGenome) : false;
    }

    public static boolean readMismatchesPastJunction(final Read read, final Junction junction, final RefGenomeInterface refGenome)
    {
        int posStart, posEnd;

        if(junction.isForward())
        {
            if(read.alignmentEnd() <= junction.Position)
                return false;

            posStart = junction.Position + 1;
            posEnd = read.alignmentEnd();
        }
        else
        {
            if(read.alignmentStart() >= junction.Position)
                return false;

            posStart = read.alignmentStart();
            posEnd = junction.Position - 1;
        }

        int refLength = posEnd - posStart + 1;

        if(refLength < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false; // insufficient for junction support anyway

        if(refLength > MIN_VARIANT_LENGTH / 2)
            return false; // should have soft-clipped

        byte[] refBases = refGenome.getBases(read.chromosome(), posStart, posEnd);

        if(refBases == null || refBases.length < refLength)
            return false;

        int readIndexStart = readIndexFromPosition(read, posStart);
        int mismatches = 0;

        for(int i = 0; i < refLength; ++i)
        {
            if(read.getBases()[readIndexStart + i] != refBases[i])
                ++mismatches;

            if(mismatches >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                return true;
        }

        return false;
    }

    public static boolean recordSoftClipsAtJunction(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() && read.alignmentEnd() == junction.Position;
        }
        else
        {
            return read.isLeftClipped() && read.alignmentStart() == junction.Position;
        }
    }

    public static int readJunctionExtensionLength(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() ? max(read.maxUnclippedEnd() - junction.Position, 0) : 0;
        }
        else
        {
            return read.isLeftClipped() ? max(junction.Position - read.minUnclippedStart(), 0) : 0;
        }
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
