package com.hartwig.hmftools.sage.vis;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.sync.FragmentData;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class ReadEvidenceRecord implements Comparable<ReadEvidenceRecord>
{
    public final SAMRecord Read;
    public final FragmentData Fragment;
    public final ReadContextMatch MatchType;
    public final QualityScores Qualities;

    private final int mVariantPosition;

    public ReadEvidenceRecord(
            final SAMRecord read, @Nullable final FragmentData fragment, final ReadContextMatch matchType,
            final QualityScores qualities, int variantPosition)
    {
        Read = read;
        Fragment = fragment;
        MatchType = matchType;
        Qualities = qualities;

        mVariantPosition = variantPosition;
    }

    private static int intCompare(int a, int b)
    {
        return a - b;
    }

    public Orientation orientation()
    {
        if(Fragment == null)
        {
            return Read.getReadNegativeStrandFlag() ? Orientation.REVERSE : Orientation.FORWARD;
        }

        boolean firstOverlaps = positionWithin(mVariantPosition, Fragment.First.getAlignmentStart(), Fragment.First.getAlignmentEnd());
        boolean secondOverlaps = positionWithin(mVariantPosition, Fragment.Second.getAlignmentStart(), Fragment.Second.getAlignmentEnd());
        Orientation firstOrientation = Fragment.First.getReadNegativeStrandFlag() ? Orientation.REVERSE : Orientation.FORWARD;
        Orientation secondOrientation = Fragment.Second.getReadNegativeStrandFlag() ? Orientation.REVERSE : Orientation.FORWARD;
        if(firstOverlaps && !secondOverlaps)
        {
            return firstOrientation;
        }

        if(!firstOverlaps && secondOverlaps)
        {
            return secondOrientation;
        }

        if(firstOverlaps && secondOverlaps)
        {
            if(firstOrientation == secondOrientation)
            {
                return firstOrientation;
            }

            return Orientation.MIXED;
        }

        if(Fragment.First.getFirstOfPairFlag())
        {
            return firstOrientation;
        }

        return secondOrientation;
    }

    @Override
    public int compareTo(final ReadEvidenceRecord other)
    {
        int orientationCompare = intCompare(orientation().SortKey, other.orientation().SortKey);
        if(orientationCompare != 0)
        {
            return orientationCompare;
        }

        if(Qualities == null && other.Qualities == null)
        {
            return -intCompare(Read.getMappingQuality(), other.Read.getMappingQuality());
        }

        if(Qualities == null)
        {
            return 1;
        }

        if(other.Qualities == null)
        {
            return -1;
        }

        if(Qualities.ModifiedQuality == other.Qualities.ModifiedQuality)
        {
            return -intCompare(Read.getMappingQuality(), other.Read.getMappingQuality());
        }

        return Qualities.ModifiedQuality > other.Qualities.ModifiedQuality ? -1 : 1;
    }

    public enum Orientation
    {
        FORWARD(-1),
        REVERSE(1),
        MIXED(0);

        public final int SortKey;

        private Orientation(int sortKey)
        {
            SortKey = sortKey;
        }
    }
}