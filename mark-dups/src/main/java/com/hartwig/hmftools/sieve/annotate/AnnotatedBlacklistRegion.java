package com.hartwig.hmftools.sieve.annotate;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class AnnotatedBlacklistRegion
{
    public static final String CSV_HEADER =
            BlacklistRegion.CSV_HEADER + "," + AnnotateStatistics.CSV_HEADER + "," + RepeatMasker.CSV_HEADER + ","
                    + AnnotateStatistics.CSV_HEADER;

    final private BlacklistRegion mBlacklistRegion;
    final private List<Pair<RepeatMasker, AnnotateStatistics>> mAnnotatedRepeatMaskers;
    final private AnnotateStatistics mStatistics;

    public AnnotatedBlacklistRegion(@NotNull final BlacklistRegion blacklistRegion)
    {
        mBlacklistRegion = blacklistRegion;
        mAnnotatedRepeatMaskers = new ArrayList<>();
        mStatistics = new AnnotateStatistics();
    }

    public void addRepeatMasker(@NotNull final RepeatMasker repeatMasker)
    {
        // TODO(m_cooper): Do this check.
        //        if (repeatMasker.getRepeatPosStart() < mBlacklistRegion.getPosStart() || repeatMasker.getRepeatPosEnd() > mBlacklistRegion.getPosEnd())
        //        {
        //            MD_LOGGER.error("Repeat masker ({}) is not contained within the blacklist region ({}).", repeatMasker.toString(), mBlacklistRegion.toString());
        //            System.exit(1);
        //        }

        mAnnotatedRepeatMaskers.add(Pair.of(repeatMasker, new AnnotateStatistics()));
    }

    @NotNull
    public BlacklistRegion getBlacklistRegion()
    {
        return mBlacklistRegion;
    }

    @NotNull
    public List<String> getCSVLines()
    {
        final String blacklistRegionFragment = mBlacklistRegion.getCSVFragment();
        final String statisticsFragment = mStatistics.getCSVFragment();
        if(mAnnotatedRepeatMaskers.isEmpty())
        {
            StringBuilder sb = new StringBuilder();
            sb.append(blacklistRegionFragment);
            sb.append(',');
            sb.append(statisticsFragment);
            sb.append(',');
            sb.append(RepeatMasker.EMPTY_CSV_FRAGMENT);
            sb.append(',');
            sb.append(AnnotateStatistics.EMPTY_CSV_FRAGMENT);
            return List.of(sb.toString());
        }

        final List<String> lines = new ArrayList<>();
        for(var annotatedRepeatMasker : mAnnotatedRepeatMaskers)
        {
            StringBuilder sb  = new StringBuilder();
            sb.append(blacklistRegionFragment);
            sb.append(',');
            sb.append(statisticsFragment);
            sb.append(',');
            sb.append(annotatedRepeatMasker.getLeft().getCSVFragment());
            sb.append(',');
            sb.append(annotatedRepeatMasker.getRight().getCSVFragment());
            lines.add(sb.toString());
        }

        return lines;
    }

    public void matchedRead(@NotNull final SAMRecord read)
    {
        mStatistics.matchedRead(read);

        // TODO(m_cooper): What if the repeat mask is outside of the blacklisted region?
        for(var annotatedRepeatMask : mAnnotatedRepeatMaskers)
        {
            RepeatMasker repeatMask = annotatedRepeatMask.getLeft();
            AnnotateStatistics repeatMaskStats = annotatedRepeatMask.getRight();
            if(read.getAlignmentStart() >= repeatMask.getRepeatPosStart() && read.getAlignmentEnd() <= repeatMask.getRepeatPosEnd())
            {
                repeatMaskStats.matchedRead(read);
            }
        }
    }
}
