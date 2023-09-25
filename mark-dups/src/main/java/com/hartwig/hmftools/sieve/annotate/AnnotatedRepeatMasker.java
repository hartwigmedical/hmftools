package com.hartwig.hmftools.sieve.annotate;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class AnnotatedRepeatMasker implements IJobRegion
{
    private final AnnotatedBlacklistRegion mParentRegion;
    private final RepeatMasker mRepeatMasker;
    private final AnnotateStatistics mStats;

    public AnnotatedRepeatMasker(@NotNull final AnnotatedBlacklistRegion parentRegion, @NotNull final RepeatMasker repeatMasker)
    {
        mParentRegion = parentRegion;
        mRepeatMasker = repeatMasker;
        mStats = new AnnotateStatistics();
    }

    @Override
    public void matchedRead(@NotNull final SAMRecord read)
    {
        mStats.matchedRead(read);
        // TODO(m_cooper): What if the repeat mask is outside of the blacklisted region?
    }

    @Override
    public ChrBaseRegion getChrBaseRegion()
    {
        return new ChrBaseRegion(
                mParentRegion.getBlacklistRegion().getChromosome(),
                mRepeatMasker.getRepeatPosStart(),
                mRepeatMasker.getRepeatPosEnd());
    }

    @NotNull
    public String getCSVFragment()
    {
        return mRepeatMasker.getCSVFragment() + ',' + mStats.getCSVFragment();
    }

    @NotNull
    public AnnotatedBlacklistRegion getParentRegion()
    {
        return mParentRegion;
    }

    @NotNull
    public RepeatMasker getRepeatMasker()
    {
        return mRepeatMasker;
    }
}
