package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.SAMRecord;

public class FragmentCoords implements Comparable<FragmentCoords>
{
    public final String ChromsomeLower;
    public final String ChromsomeUpper;
    public final int PositionLower;
    public final int PositionUpper;
    public final Orientation OrientLower;
    public final Orientation OrientUpper;
    public final boolean ReadIsLower; // the read which created these coordinates
    public final boolean SupplementarySourced;
    public final boolean UnmappedSourced;
    public final boolean IsForward; // forward = F1R2, reverse is F2R1 - relates to collapsing and dual-strand classification
    public final String Key;

    public FragmentCoords(
            final String chromsomeLower, final String chromsomeUpper, final int positionLower, final int positionUpper,
            final Orientation orientLower, final Orientation orientUpper,
            boolean readIsLower, boolean isForwardOrientation, boolean supplementarySourced, final boolean unmappedSourced)
    {
        ChromsomeLower = chromsomeLower;
        ChromsomeUpper = chromsomeUpper;
        PositionLower = positionLower;
        PositionUpper = positionUpper;
        OrientLower = orientLower;
        OrientUpper = orientUpper;
        ReadIsLower = readIsLower;
        IsForward = isForwardOrientation;
        SupplementarySourced = supplementarySourced;
        UnmappedSourced = unmappedSourced;

        String coordinateLower = formCoordinate(ChromsomeLower, positionLower, orientLower.isForward());

        if(positionUpper == NO_POSITION)
        {
            if(unmappedSourced)
                Key = format("%s_U", coordinateLower);
            else
                Key = format("%s", coordinateLower);
        }
        else
        {
            String coordinateUpper = formCoordinate(ChromsomeUpper, positionUpper, orientUpper.isForward());
            Key = format("%s_%s", coordinateLower, coordinateUpper);
        }
    }

    public boolean isSingleRead() { return PositionUpper == NO_POSITION; }

    public static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }

    @Override
    public int compareTo(final FragmentCoords other)
    {
        return Key.compareTo(other.Key);
    }

    @Override
    public int hashCode() { return Key.hashCode(); }

    public String toString() { return Key; }

    public static FragmentCoords fromRead(final SAMRecord read)
    {
        Orientation readOrient;
        String readChromosome;
        int readPosition;

        Orientation mateOrient = null;
        String mateChromosome = NO_CHROMOSOME_NAME;
        int matePosition = NO_POSITION;

        boolean isPaired = read.getReadPairedFlag();
        boolean isSupplementary = read.getSupplementaryAlignmentFlag();

        boolean isUnmapped = read.getReadUnmappedFlag();

        if(isPaired && !read.getMateUnmappedFlag())
        {
            mateChromosome = read.getMateReferenceName();
            mateOrient = read.getMateNegativeStrandFlag() ? Orientation.REVERSE : Orientation.FORWARD;

            String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            matePosition = getFivePrimeUnclippedPosition(read.getMateAlignmentStart(), mateCigar, mateOrient.isForward());
        }

        if(isSupplementary)
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            readChromosome = suppData.Chromosome;
            readOrient = Orientation.fromChar(suppData.Strand);
            readPosition = getFivePrimeUnclippedPosition(suppData.Position, suppData.Cigar, readOrient.isForward());
        }
        else if(!isPaired || !isUnmapped)
        {
            readChromosome = read.getReferenceName();
            readOrient = read.getReadNegativeStrandFlag() ? Orientation.REVERSE : Orientation.FORWARD;

            readPosition = read.getCigar() != null ?
                    getFivePrimeUnclippedPosition(read) :
                    getFivePrimeUnclippedPosition(read.getAlignmentStart(), read.getCigarString(), readOrient.isForward());
        }
        else
        {
            // use mate's coords for unmapped reads
            readChromosome = mateChromosome;
            readOrient = mateOrient;
            readPosition = matePosition;
        }

        boolean readIsLower;

        if(!isPaired || isUnmapped || read.getMateUnmappedFlag())
        {
            readIsLower = isUnmapped;

            return new FragmentCoords(
                    readChromosome, NO_CHROMOSOME_NAME, readPosition, NO_POSITION, readOrient, readOrient,
                    readIsLower, true, isSupplementary, isUnmapped);
        }

        if(readChromosome.equals(mateChromosome))
        {
            readIsLower = readPosition <= matePosition;
        }
        else
        {
            readIsLower = read.getReferenceIndex() < read.getMateReferenceIndex();
        }

        boolean isForwardOrientation = readIsLower == read.getFirstOfPairFlag();

        if(readIsLower)
        {
            return new FragmentCoords(
                    readChromosome, mateChromosome, readPosition, matePosition, readOrient, mateOrient,
                    readIsLower, isForwardOrientation, isSupplementary, isUnmapped);
        }
        else
        {
            return new FragmentCoords(
                    mateChromosome, readChromosome, matePosition, readPosition, mateOrient, readOrient,
                    readIsLower, isForwardOrientation, isSupplementary, isUnmapped);
        }
    }
}
