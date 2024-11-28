package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;

import org.checkerframework.checker.units.qual.K;

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
    public final Orientation FragmentOrient; // forward = F1R2, reverse is F2R1 - relates to collapsing and dual-strand classification
    public final String Key;
    public final String KeyOriented; // includes fragment orientation

    public FragmentCoords(
            final String chromsomeLower, final String chromsomeUpper, final int positionLower, final int positionUpper,
            final Orientation orientLower, final Orientation orientUpper, final Orientation fragmentOrientation,
            boolean readIsLower, boolean isSupplementary, boolean isUnmapped, boolean keyByFragmentOrientation)
    {
        ChromsomeLower = chromsomeLower;
        ChromsomeUpper = chromsomeUpper;
        PositionLower = positionLower;
        PositionUpper = positionUpper;
        OrientLower = orientLower;
        OrientUpper = orientUpper;
        FragmentOrient = fragmentOrientation;
        ReadIsLower = readIsLower;
        SupplementarySourced = isSupplementary;
        UnmappedSourced = isUnmapped;

        String coordinateLower = formCoordinate(ChromsomeLower, positionLower, orientLower.isForward());

        if(positionUpper == NO_POSITION)
        {
            if(isUnmapped)
                Key = format("%s_U", coordinateLower);
            else
                Key = format("%s", coordinateLower);

            KeyOriented = Key;
        }
        else
        {
            String coordinateUpper = formCoordinate(ChromsomeUpper, positionUpper, orientUpper.isForward());

            Key = format("%s_%s", coordinateLower, coordinateUpper);
            KeyOriented = (keyByFragmentOrientation && FragmentOrient.isReverse()) ? format("%s_N", Key) : Key;
        }
    }

    public boolean isSingleRead() { return PositionUpper == NO_POSITION; }
    public boolean forwardFragment() { return FragmentOrient.isForward(); }

    public static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }

    public static FragmentCoords fromRead(final SAMRecord read) { return fromRead(read, false); }

    public static FragmentCoords fromRead(final SAMRecord read, boolean useFragmentOrientation)
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
            mateOrient = read.getMateNegativeStrandFlag() ? REVERSE : FORWARD;

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
            readOrient = read.getReadNegativeStrandFlag() ? REVERSE : FORWARD;

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
            readIsLower = !isUnmapped;

            return new FragmentCoords(
                    readChromosome, NO_CHROMOSOME_NAME, readPosition, NO_POSITION, readOrient, readOrient, FORWARD,
                    readIsLower, isSupplementary, isUnmapped, false);
        }

        if(readChromosome.equals(mateChromosome))
        {
            readIsLower = readPosition <= matePosition;
        }
        else
        {
            readIsLower = read.getReferenceIndex() < read.getMateReferenceIndex();
        }

        Orientation fragmentOrientation = (readIsLower == read.getFirstOfPairFlag()) ? FORWARD : REVERSE;

        if(readIsLower)
        {
            return new FragmentCoords(
                    readChromosome, mateChromosome, readPosition, matePosition, readOrient, mateOrient, fragmentOrientation,
                    readIsLower, isSupplementary, isUnmapped, useFragmentOrientation);
        }
        else
        {
            return new FragmentCoords(
                    mateChromosome, readChromosome, matePosition, readPosition, mateOrient, readOrient, fragmentOrientation,
                    readIsLower, isSupplementary, isUnmapped, useFragmentOrientation);
        }
    }

    @Override
    public int compareTo(final FragmentCoords other)
    {
        return KeyOriented.compareTo(other.KeyOriented);
    }

    @Override
    public boolean equals(@Nullable Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof FragmentCoords))
            return false;

        FragmentCoords fragCoords = (FragmentCoords)other;
        return KeyOriented.compareTo(fragCoords.KeyOriented) == 0;
    }

    @Override
    public int hashCode() { return KeyOriented.hashCode(); }

    public String toString() { return KeyOriented; }
}
