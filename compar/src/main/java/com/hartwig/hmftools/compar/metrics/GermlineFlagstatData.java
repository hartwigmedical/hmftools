package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineFlagstatData implements ComparableItem
{
    public final BamFlagStats mFlagstat;

    public GermlineFlagstatData(final BamFlagStats flagstat)
    {
        mFlagstat = flagstat;
    }

    @Override
    public Category category()
    {
        return GERMLINE_FLAGSTAT;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", mFlagstat.mappedProportion()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final GermlineFlagstatData otherData = (GermlineFlagstatData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_MAPPED_PROPORTION, mFlagstat.mappedProportion(), otherData.mFlagstat.mappedProportion(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, includeMatches);
    }
}
