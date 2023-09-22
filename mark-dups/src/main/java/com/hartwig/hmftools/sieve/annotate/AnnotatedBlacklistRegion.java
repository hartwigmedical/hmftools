package com.hartwig.hmftools.sieve.annotate;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class AnnotatedBlacklistRegion
{
    public static final String CSV_HEADER = BlacklistRegion.CSV_HEADER + "," + RepeatMasker.CSV_HEADER;

    final BlacklistRegion mBlacklistRegion;
    final List<RepeatMasker> mRepeatMaskers;

    public AnnotatedBlacklistRegion(@NotNull final BlacklistRegion blacklistRegion)
    {
        mBlacklistRegion = blacklistRegion;
        mRepeatMaskers = new ArrayList<>();
    }

    public void addRepeatMasker(@NotNull final RepeatMasker repeatMasker)
    {
        // TODO(m_cooper): Do this check.
        //        if (repeatMasker.getRepeatPosStart() < mBlacklistRegion.getPosStart() || repeatMasker.getRepeatPosEnd() > mBlacklistRegion.getPosEnd())
        //        {
        //            MD_LOGGER.error("Repeat masker ({}) is not contained within the blacklist region ({}).", repeatMasker.toString(), mBlacklistRegion.toString());
        //            System.exit(1);
        //        }

        mRepeatMaskers.add(repeatMasker);
    }

    @NotNull
    public BlacklistRegion getBlacklistRegion()
    {
        return mBlacklistRegion;
    }

    @NotNull
    public List<RepeatMasker> getRepeatMaskers()
    {
        return mRepeatMaskers;
    }

    @NotNull
    public List<String> getCSVLines()
    {
        final String blacklistRegionFragment = mBlacklistRegion.getCSVFragment();
        if(mRepeatMaskers.isEmpty())
        {
            final StringBuilder sb = new StringBuilder();
            sb.append(blacklistRegionFragment);
            sb.append(',');
            sb.append(RepeatMasker.EMPTY_CSV_FRAGMENT);
            return List.of(sb.toString());
        }

        final List<String> lines = new ArrayList<>();
        for(RepeatMasker repeatMasker : mRepeatMaskers)
        {
            final StringBuilder sb = new StringBuilder();
            sb.append(blacklistRegionFragment);
            sb.append(',');
            sb.append(repeatMasker.getCSVFragment());
            lines.add(sb.toString());
        }

        return lines;
    }
}
