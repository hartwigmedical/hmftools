package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.loadBamMetricsSummary;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class BamMetricsComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public final CategoryType mCategory;
    private final List<Integer> mComparisonPercentages;

    protected static final List<Integer> TUMOR_COVERAGE_PERCENTAGES = List.of(30, 60);
    protected static final List<Integer> GERMLINE_COVERAGE_PERCENTAGES = List.of(10, 20);

    public BamMetricsComparer(final CategoryType categoryType, final ComparConfig config)
    {
        mCategory = categoryType;
        mComparisonPercentages = categoryType == TUMOR_BAM_METRICS ? TUMOR_COVERAGE_PERCENTAGES : GERMLINE_COVERAGE_PERCENTAGES;
        mConfig = config;
    }

    @Override
    public CategoryType category()
    {
        return mCategory;
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

        for(Integer coverage : mComparisonPercentages)
        {
            thresholds.addFieldThreshold(coverageString(coverage), 0.03, 0);
        }
    }

    protected static String coverageString(final int coverage) { return format("Percentage%dX", coverage); }

    @Override
    public List<String> comparedFieldNames()
    {
        List<String> fieldNames = Lists.newArrayList(FLD_DUPLICATE_PERCENTAGE);
        for(Integer coverage : mComparisonPercentages)
        {
            fieldNames.add(coverageString(coverage));
        }

        return fieldNames;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            BamMetricSummary metrics = loadBamMetricsSummary(sampleId, fileSources.TumorBamMetrics);
            comparableItems.add(new BamMetricsData(mCategory, metrics, mComparisonPercentages));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load tumor BAM metrics data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }
}
