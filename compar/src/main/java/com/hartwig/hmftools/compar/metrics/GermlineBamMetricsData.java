package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineBamMetricsData implements ComparableItem
{
    public final BamMetricsSummary Metrics;

    protected static final String FLD_PERCENTAGE_10X = "Percentage10X";
    protected static final String FLD_PERCENTAGE_20X = "Percentage20X";

    public GermlineBamMetricsData(final BamMetricsSummary metrics)
    {
        Metrics = metrics;
    }

    @Override
    public Category category()
    {
        return GERMLINE_BAM_METRICS;
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
        values.add(format("%.2f", Metrics.duplicatePercent()));
        values.add(format("%.2f", Metrics.coveragePercent(10)));
        values.add(format("%.2f", Metrics.coveragePercent(20)));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean isPass() {
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
        final GermlineBamMetricsData otherData = (GermlineBamMetricsData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_DUPLICATE_PERCENTAGE, Metrics.duplicatePercent(), otherData.Metrics.duplicatePercent(), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_10X, Metrics.coveragePercent(10), otherData.Metrics.coveragePercent(10), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_20X, Metrics.coveragePercent(20), otherData.Metrics.coveragePercent(20), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
