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

import htsjdk.samtools.SAMRecord;

public class FragmentCoords implements Comparable<FragmentCoords>
{
    public final String ChromsomeLower;
    public final String ChromsomeUpper;
    public final int PositionLower;
    public final int PositionUpper;
    public final Orientation OrientLower;
    public final Orientation OrientUpper;
    public final Orientation FragmentOrient; // forward = F1R2, reverse is F2R1 - relates to collapsing and dual-strand classification

    // info about the read used to create these coordinates
    public final boolean ReadIsLower; // read set the power fragment position
    public final SupplementaryReadInfo SuppReadInfo; // set if sourced from a supplementary read
    public final boolean UnmappedSourced;

    public final String Key;
    public final String KeyOriented; // includes fragment orientation

    public FragmentCoords(
            final String chromsomeLower, final String chromsomeUpper, final int positionLower, final int positionUpper,
            final Orientation orientLower, final Orientation orientUpper, final Orientation fragmentOrientation,
            boolean readIsLower, final SupplementaryReadInfo suppReadInfo, boolean isUnmapped, boolean keyByFragmentOrientation)
    {
        ChromsomeLower = chromsomeLower;
        ChromsomeUpper = chromsomeUpper;
        PositionLower = positionLower;
        PositionUpper = positionUpper;
        OrientLower = orientLower;
        OrientUpper = orientUpper;
        FragmentOrient = fragmentOrientation;

        ReadIsLower = readIsLower;
        SuppReadInfo = suppReadInfo;
        UnmappedSourced = isUnmapped;

        String coordinateLower = formCoordinate(ChromsomeLower, positionLower, orientLower.isForward());

        if(positionUpper == NO_POSITION)
        {
            Key = isUnmapped ? format("%s_U", coordinateLower) : format("%s", coordinateLower);
            KeyOriented = Key;
        }
        else
        {
            char lowerUpperChar = ReadIsLower ? 'L' : 'U';
            String coordinateUpper = formCoordinate(ChromsomeUpper, positionUpper, orientUpper.isForward());

            Key = format("%c_%s_%s", lowerUpperChar, coordinateLower, coordinateUpper);
            KeyOriented = (keyByFragmentOrientation && FragmentOrient.isReverse()) ? format("%s_N", Key) : Key;
        }
    }

    public boolean isSingleRead() { return PositionUpper == NO_POSITION; }
    public boolean forwardFragment() { return FragmentOrient.isForward(); }

    public int readPosition()
    {
        if(SuppReadInfo != null)
            return SuppReadInfo.UnclippedPosition;
        else
            return ReadIsLower || UnmappedSourced ? PositionLower : PositionUpper;
    }

    public Orientation readOrientation()
    {
        if(SuppReadInfo != null)
            return SuppReadInfo.Orient;
        else
            return ReadIsLower || UnmappedSourced ? OrientLower : OrientUpper;
    }

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

        SupplementaryReadInfo suppReadInfo = null;

        if(!isPaired || !isUnmapped)
        {
            readChromosome = read.getReferenceName();
            readOrient = read.getReadNegativeStrandFlag() ? REVERSE : FORWARD;

            readPosition = read.getCigar() != null ?
                    getFivePrimeUnclippedPosition(read) :
                    getFivePrimeUnclippedPosition(read.getAlignmentStart(), read.getCigarString(), readOrient.isForward());

            if(isSupplementary)
            {
                // use the primaries coordinates for the fragment but retain the read's info for use in the read cache
                suppReadInfo = new SupplementaryReadInfo(readPosition, readOrient);

                SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

                readChromosome = suppData.Chromosome;
                readOrient = Orientation.fromChar(suppData.Strand);
                readPosition = getFivePrimeUnclippedPosition(suppData.Position, suppData.Cigar, readOrient.isForward());
            }
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
                    readIsLower, suppReadInfo, isUnmapped, false);
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
                    readIsLower, suppReadInfo, isUnmapped, useFragmentOrientation);
        }
        else
        {
            return new FragmentCoords(
                    mateChromosome, readChromosome, matePosition, readPosition, mateOrient, readOrient, fragmentOrientation,
                    readIsLower, suppReadInfo, isUnmapped, useFragmentOrientation);
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
