package com.hartwig.hmftools.sieve.annotate;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class AnnotatedBlacklistRegion implements IJobRegion
{
    public static final String CSV_HEADER =
            BlacklistRegion.CSV_HEADER + "," + AnnotateStatistics.CSV_HEADER + "," + RepeatMasker.CSV_HEADER + ","
                    + AnnotateStatistics.CSV_HEADER;

    final private BlacklistRegion mBlacklistRegion;
    final private AnnotateStatistics mStats;

    final private List<AnnotatedRepeatMasker> mAnnotatedRepeatMaskers;

    public AnnotatedBlacklistRegion(@NotNull final BlacklistRegion blacklistRegion)
    {
        mBlacklistRegion = blacklistRegion;
        mAnnotatedRepeatMaskers = new ArrayList<>();
        mStats = new AnnotateStatistics();
    }

    public void addRepeatMasker(@NotNull final RepeatMasker repeatMasker)
    {
        // TODO(m_cooper): Do this check.
        //        if (repeatMasker.getRepeatPosStart() < mBlacklistRegion.getPosStart() || repeatMasker.getRepeatPosEnd() > mBlacklistRegion.getPosEnd())
        //        {
        //            MD_LOGGER.error("Repeat masker ({}) is not contained within the blacklist region ({}).", repeatMasker.toString(), mBlacklistRegion.toString());
        //            System.exit(1);
        //        }

        mAnnotatedRepeatMaskers.add(new AnnotatedRepeatMasker(this, repeatMasker));
    }

    @NotNull
    public BlacklistRegion getBlacklistRegion()
    {
        return mBlacklistRegion;
    }

    public List<IJobRegion> getJobRegions()
    {
        List<IJobRegion> output = new ArrayList<>();
        output.add(this);
        output.addAll(mAnnotatedRepeatMaskers);
        return output;
    }

    @NotNull
    public List<String> getCSVLines()
    {
        final String blacklistRegionFragment = mBlacklistRegion.getCSVFragment();
        final String statsFragment = mStats.getCSVFragment();
        if(mAnnotatedRepeatMaskers.isEmpty())
        {
            final String line = blacklistRegionFragment + ',' + statsFragment + ',' + RepeatMasker.EMPTY_CSV_FRAGMENT + ','
                    + AnnotateStatistics.EMPTY_CSV_FRAGMENT;
            return List.of(line);
        }

        final List<String> lines = new ArrayList<>();
        for(var annotatedRepeatMasker : mAnnotatedRepeatMaskers)
        {
            final String line = blacklistRegionFragment + ',' + statsFragment + ',' + annotatedRepeatMasker.getCSVFragment();

            lines.add(line);
        }

        return lines;
    }

    @Override
    public ChrBaseRegion getChrBaseRegion()
    {
        return new ChrBaseRegion(mBlacklistRegion.getChromosome(), mBlacklistRegion.getPosStart(), mBlacklistRegion.getPosEnd());
    }

    @Override
    public void matchedRead(@NotNull final SAMRecord read)
    {
        mStats.matchedRead(read);
    }
}
