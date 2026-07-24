package com.hartwig.hmftools.compar.metrics;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class BamMetricsData implements ComparableItem
{
    private final CategoryType mCategory;
    private final BamMetricSummary mMetrics;

    public BamMetricsData(
            final CategoryType category, final BamMetricSummary metrics)
    {
        mMetrics = metrics;
        mCategory = category;
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
