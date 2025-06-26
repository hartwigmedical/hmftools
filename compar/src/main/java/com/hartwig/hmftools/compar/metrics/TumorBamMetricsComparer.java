package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.loadBamMetricsSummary;
import static com.hartwig.hmftools.compar.metrics.TumorBamMetricsData.FLD_PERCENTAGE_30X;
import static com.hartwig.hmftools.compar.metrics.TumorBamMetricsData.FLD_PERCENTAGE_60X;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SampleFileSources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class TumorBamMetricsComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public TumorBamMetricsComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return TUMOR_BAM_METRICS;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_DUPLICATE_PERCENTAGE, DUPLICATE_PERCENTAGE_ABS_THRESHOLD, DUPLICATE_PERCENTAGE_PCT_THRESHOLD);
        thresholds.addFieldThreshold(FLD_PERCENTAGE_30X, 0.03, 0);
        thresholds.addFieldThreshold(FLD_PERCENTAGE_60X, 0.03, 0);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_DUPLICATE_PERCENTAGE, FLD_PERCENTAGE_30X, FLD_PERCENTAGE_60X);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final SampleFileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            BamMetricsSummary metrics = loadBamMetricsSummary(sampleId, fileSources.tumorBamMetrics());
            comparableItems.add(new TumorBamMetricsData(metrics));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load tumor BAM metrics data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }
}
