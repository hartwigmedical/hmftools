package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getThreePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.SAMRecord;

public class FragmentCoords implements Comparable<FragmentCoords>
{
    // represents all the fragment property of a read or pair of reads
    // primary purpose is to find duplicates, by a match in their Key, which is defined as:
    // A_B_C_D where
    // A: LowerCoordinate, where coordinate is chrmosome:unclipped_position and :R if reverse orientation
    // B: UpperCoordinate, or the read's unclipped end if unpaired
    // C: whether the read which builds the coordinates matched the lower or upper position if a paired read
    // D: if the read is unmapped or supplementary

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
    public final boolean Unpaired;

    public final String Key; // includes fragment orientation if applicable

    private final String mKeyNonOriented; // only used for UMI collapsing, and set as required

    public FragmentCoords(
            final String chromsomeLower, final String chromsomeUpper, final int positionLower, final int positionUpper,
            final Orientation orientLower, final Orientation orientUpper, final Orientation fragmentOrientation, boolean readIsLower,
            final SupplementaryReadInfo suppReadInfo, boolean isUnmapped, boolean keyByFragmentOrientation, boolean unpaired)
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
        Unpaired = unpaired;

        if(Unpaired)
        {
            String coordinate = format("%s:%d_%d", ChromsomeLower, PositionLower, PositionUpper);

            if(ReadIsLower ? OrientLower.isReverse() : OrientUpper.isReverse())
                coordinate = format("%s_R", coordinate);

            Key = suppReadInfo != null ? format("%s_S", coordinate) : coordinate;
            mKeyNonOriented = null;
            return;
        }

        String coordinateLower = formCoordinate(ChromsomeLower, positionLower, orientLower.isForward());

        if(positionUpper == NO_POSITION)
        {
            if(isUnmapped)
                Key = format("%s_U", coordinateLower);
            else if(suppReadInfo != null)
                Key = format("%s_S", coordinateLower);
            else
                Key = coordinateLower;

            mKeyNonOriented = null;
        }
        else
        {
            String coordinateUpper = formCoordinate(ChromsomeUpper, positionUpper, orientUpper.isForward());

            String readInfo = ReadIsLower ? "L" : "U";
            if(SuppReadInfo != null)
                readInfo += "_S";

            String keyNonOriented = format("%s_%s_%s", coordinateLower, coordinateUpper, readInfo);
            if(keyByFragmentOrientation && FragmentOrient.isReverse())
            {
                Key = format("%s_N", keyNonOriented);
                mKeyNonOriented = keyNonOriented;
            }
            else
            {
                Key = keyNonOriented;
                mKeyNonOriented = null;
            }
        }
    }

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

    public String keyNonOriented() { return mKeyNonOriented != null ? mKeyNonOriented : Key; }

    private static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s:%d", chromosome, position) : format("%s:%d:R", chromosome, position);
    }

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

            // if mate CIGAR isn't present, continue assuming a fully-aligned mate read - only supplementaries conditionally use this logic
            if(mateCigar == null)
                mateCigar = format("%dM", read.getReadBases().length);

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

                if(!isPaired)
                {
                    mateChromosome = readChromosome;
                    mateOrient = readOrient.opposite();
                    matePosition = getThreePrimeUnclippedPosition(suppData.Position, suppData.Cigar, readOrient.isForward());
                }
            }
            else if(!isPaired)
            {
                mateChromosome = readChromosome;
                mateOrient = readOrient.opposite();

                matePosition = read.getCigar() != null ?
                        getThreePrimeUnclippedPosition(read) :
                        getThreePrimeUnclippedPosition(read.getAlignmentStart(), read.getCigarString(), readOrient.isForward());
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

        if(isUnmapped || (isPaired && read.getMateUnmappedFlag()))
        {
            readIsLower = !isUnmapped;

            return new FragmentCoords(
                    readChromosome, NO_CHROMOSOME_NAME, readPosition, NO_POSITION, readOrient, readOrient, FORWARD,
                    readIsLower, suppReadInfo, isUnmapped, false, !isPaired);
        }

        // the following determination of which of the primary reads is consider 'lower' in the fragment must give the same result
        // from each of their perspectives
        if(!readChromosome.equals(mateChromosome))
        {
            readIsLower = read.getReferenceIndex() < read.getMateReferenceIndex();
        }
        else if(readPosition != matePosition)
        {
            readIsLower = readPosition < matePosition;
        }
        else
        {
            readIsLower = read.getFirstOfPairFlag();
        }

        Orientation fragmentOrientation = (!isPaired || (readIsLower == read.getFirstOfPairFlag())) ? FORWARD : REVERSE;

        if(readIsLower)
        {
            return new FragmentCoords(
                    readChromosome, mateChromosome, readPosition, matePosition, readOrient, mateOrient, fragmentOrientation,
                    readIsLower, suppReadInfo, isUnmapped, useFragmentOrientation, !isPaired);
        }
        else
        {
            return new FragmentCoords(
                    mateChromosome, readChromosome, matePosition, readPosition, mateOrient, readOrient, fragmentOrientation,
                    readIsLower, suppReadInfo, isUnmapped, useFragmentOrientation, !isPaired);
        }
    }

    public FragmentCoords withFragmentOrientation(final Orientation fragmentOrientation)
    {
        if(FragmentOrient == fragmentOrientation)
            return this;

        return new FragmentCoords(ChromsomeLower, ChromsomeUpper, PositionLower, PositionUpper, OrientLower, OrientUpper, fragmentOrientation, ReadIsLower, SuppReadInfo, UnmappedSourced, true, UnmappedSourced);
    }

    @Override
    public int compareTo(final FragmentCoords other)
    {
        return Key.compareTo(other.Key);
    }

    @Override
    public boolean equals(@Nullable Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof FragmentCoords))
            return false;

        FragmentCoords fragCoords = (FragmentCoords)other;
        return Key.compareTo(fragCoords.Key) == 0;
    }

    @Override
    public int hashCode() { return Key.hashCode(); }

    public String toString() { return Key; }
}
