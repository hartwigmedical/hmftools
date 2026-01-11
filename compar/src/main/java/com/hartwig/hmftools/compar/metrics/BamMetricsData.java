package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class BamMetricsData implements ComparableItem
{
    private final CategoryType mCategory;
    private final BamMetricSummary mMetrics;
    private final List<Integer> mComparisonPercentages;

    public BamMetricsData(
            final CategoryType category, final BamMetricSummary metrics, final List<Integer> comparisonPercentages)
    {
        mMetrics = metrics;
        mCategory = category;
        mComparisonPercentages = comparisonPercentages;
    }

    public BamMetricSummary metrics() { return mMetrics; }

    @Override
    public CategoryType category() { return mCategory; }

    @Override
    public String key() { return ""; }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", mMetrics.duplicatePercent()));

        for(Integer coverage : mComparisonPercentages)
        {
            values.add(format("%.2f", mMetrics.coveragePercent(coverage)));
        }

        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean isPass() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches)
    {
        final BamMetricsData otherData = (BamMetricsData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_DUPLICATE_PERCENTAGE, mMetrics.duplicatePercent(), otherData.mMetrics.duplicatePercent(), thresholds);

        for(Integer coverage : mComparisonPercentages)
        {
            String coverageStr = format("Percentage%dX", coverage);
            checkDiff(diffs, coverageStr, mMetrics.coveragePercent(coverage), otherData.mMetrics.coveragePercent(coverage), thresholds);
        }

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
