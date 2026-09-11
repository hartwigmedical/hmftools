package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.metrics.BamMetricsComparer.coverageString;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class BamMetricsData extends ComparableItem
{
    private final CategoryType mCategory;
    private final BamMetricSummary mMetrics;

    public BamMetricsData(
            final CategoryType category, final BamMetricSummary metrics,
            final List<Integer> comparisonPercentages, final List<FieldInfo> fields)
    {
        mMetrics = metrics;
        mCategory = category;

        addDoubleValue(FLD_DUPLICATE_PERCENTAGE, metrics.duplicatePercent(), fields);

        for(Integer coverage : comparisonPercentages)
        {
            String coverageStr = coverageString(coverage);
            double coveragePerc = metrics().coveragePercent(coverage);
            addDoubleValue(coverageStr, coveragePerc, fields);
        }
    }

    public BamMetricSummary metrics() { return mMetrics; }

    @Override
    public CategoryType category() { return mCategory; }

    @Override
    public String key() { return ""; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
