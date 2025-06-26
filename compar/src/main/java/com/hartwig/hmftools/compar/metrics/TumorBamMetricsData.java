package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.TUMOR_BAM_METRICS;
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

public class TumorBamMetricsData implements ComparableItem
{
    public final BamMetricsSummary Metrics;

    protected static final String FLD_PERCENTAGE_30X = "Percentage30X";
    protected static final String FLD_PERCENTAGE_60X = "Percentage60X";

    public TumorBamMetricsData(final BamMetricsSummary metrics)
    {
        Metrics = metrics;
    }

    @Override
    public Category category()
    {
        return TUMOR_BAM_METRICS;
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
        values.add(format("%.2f", Metrics.coveragePercent(30)));
        values.add(format("%.2f", Metrics.coveragePercent(60)));
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
        final TumorBamMetricsData otherData = (TumorBamMetricsData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_DUPLICATE_PERCENTAGE, Metrics.duplicatePercent(), otherData.Metrics.duplicatePercent(), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_30X, Metrics.coveragePercent(30), otherData.Metrics.coveragePercent(30), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_60X, Metrics.coveragePercent(60), otherData.Metrics.coveragePercent(60), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
