package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class TumorBamMetricsData implements ComparableItem
{
    public final WGSMetrics Metrics;

    protected static final String FLD_PERCENTAGE_30X = "Percentage30X";
    protected static final String FLD_PERCENTAGE_60X = "Percentage60X";

    public TumorBamMetricsData(final WGSMetrics metrics)
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
        values.add(format("%.2f", Metrics.pctExcDupe()));
        values.add(format("%.2f", Metrics.coverage30xPercentage()));
        values.add(format("%.2f", Metrics.coverage60xPercentage()));
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
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final TumorBamMetricsData otherData = (TumorBamMetricsData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_DUPLICATE_PERCENTAGE, Metrics.pctExcDupe(), otherData.Metrics.pctExcDupe(), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_30X, Metrics.coverage30xPercentage(), otherData.Metrics.coverage30xPercentage(), thresholds);
        checkDiff(diffs, FLD_PERCENTAGE_60X, Metrics.coverage60xPercentage(), otherData.Metrics.coverage60xPercentage(), thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
