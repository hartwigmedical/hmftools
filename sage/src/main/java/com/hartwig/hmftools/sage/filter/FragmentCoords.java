package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.getThreePrimeUnclippedPosition;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;

import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class FragmentCoords
{
    private final int mMaxCounts;
    private final Set<String> mLowerCoords;
    private final Set<String> mUpperCoords;

    public FragmentCoords(int maxCounts)
    {
        mMaxCounts = maxCounts;
        mLowerCoords = Sets.newHashSetWithExpectedSize(maxCounts);
        mUpperCoords = Sets.newHashSetWithExpectedSize(maxCounts);
    }

    public boolean atCapacity() { return mLowerCoords.size() >= mMaxCounts && mUpperCoords.size() >= mMaxCounts; }

    public void addRead(final SAMRecord read, @Nullable final SAMRecord mateRead)
    {
        if(atCapacity())
            return;

        int fivePrimeUnclippedPosition = getFivePrimeUnclippedPosition(read);
        int threePrimeUnclippedPosition = getThreePrimeUnclippedPosition(read);

        String fivePrimeReadCoord = formLocalCoordinate(fivePrimeUnclippedPosition, !read.getReadNegativeStrandFlag());
        String threePrimeReadCoord = formLocalCoordinate(threePrimeUnclippedPosition, !read.getReadNegativeStrandFlag());

        if(!read.getReadPairedFlag() || read.getMateUnmappedFlag())
        {

            mLowerCoords.add(fivePrimeUnclippedPosition < threePrimeUnclippedPosition ? fivePrimeReadCoord : threePrimeReadCoord);
            mUpperCoords.add(fivePrimeUnclippedPosition < threePrimeUnclippedPosition ? threePrimeReadCoord : fivePrimeReadCoord);
            return;
        }

        int mateUnclippedPosition;
        boolean mateForwardStrand;

        if(mateRead != null)
        {
            mateForwardStrand = !mateRead.getReadNegativeStrandFlag();
            mateUnclippedPosition = getFivePrimeUnclippedPosition(mateRead);
        }
        else
        {
            mateForwardStrand = !read.getMateNegativeStrandFlag();
            String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

            if(mateCigar == null)
            {
                mLowerCoords.add(fivePrimeReadCoord);
                return;
            }

            mateUnclippedPosition = getFivePrimeUnclippedPosition(read.getMateAlignmentStart(), mateCigar, mateForwardStrand);
        }

        boolean sameChromosome = read.getReferenceIndex() == read.getMateReferenceIndex();

        String mateCoord = sameChromosome ?
                formLocalCoordinate(mateUnclippedPosition, mateForwardStrand) :
                formRemoteCoordinate(read.getMateReferenceName(), mateUnclippedPosition, mateForwardStrand);

        if(sameChromosome && fivePrimeUnclippedPosition <= mateUnclippedPosition
        || !sameChromosome && read.getReferenceIndex() < read.getMateReferenceIndex())
        {
            mLowerCoords.add(fivePrimeReadCoord);
            mUpperCoords.add(mateCoord);
        }
        else
        {
            mLowerCoords.add(mateCoord);
            mUpperCoords.add(fivePrimeReadCoord);
        }
    }

    private static String formLocalCoordinate(final int position, final boolean isForward)
    {
        return isForward ? String.valueOf(position) : format("%d_R", position);
    }

    private static String formRemoteCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }

    public String toString() { return format("counts(lower=%d upper=%d)", mLowerCoords.size(), mUpperCoords.size()); }

    @VisibleForTesting
    public int lowerCount() { return mLowerCoords.size(); }
    public int upperCount() { return mUpperCoords.size(); }
}
