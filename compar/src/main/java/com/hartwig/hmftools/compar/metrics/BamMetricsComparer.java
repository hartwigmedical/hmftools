package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.TUMOR_BAM_METRICS;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.loadBamMetricsSummary;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class BamMetricsComparer extends ItemComparer
{
    protected final List<String> mFieldNames;

    public final CategoryType mCategory;
    private final List<Integer> mComparisonPercentages;

    protected static final List<Integer> TUMOR_COVERAGE_PERCENTAGES = List.of(30, 60);
    protected static final List<Integer> GERMLINE_COVERAGE_PERCENTAGES = List.of(10, 20);

    public BamMetricsComparer(final CategoryType categoryType, final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mCategory = categoryType;
        mComparisonPercentages = categoryType == TUMOR_BAM_METRICS ? TUMOR_COVERAGE_PERCENTAGES : GERMLINE_COVERAGE_PERCENTAGES;

        mFieldNames = Lists.newArrayList();
        mFieldNames.add(FLD_DUPLICATE_PERCENTAGE);

        mFields.add(new FieldInfo(
                FLD_DUPLICATE_PERCENTAGE,
                getOrMakeFieldCheck(
                        fieldCheckMap, FLD_DUPLICATE_PERCENTAGE, DUPLICATE_PERCENTAGE_ABS_THRESHOLD, DUPLICATE_PERCENTAGE_PCT_THRESHOLD),
                "%.3f"));

        for(Integer coverage : mComparisonPercentages)
        {
            String coverageStr = coverageString(coverage);
            mFieldNames.add(coverageStr);

            mFields.add(new FieldInfo(
                    coverageStr,
                    getOrMakeFieldCheck(fieldCheckMap, coverageStr, 0.03, null),
                    "%.3f"));
        }
    }

    @Override
    public CategoryType category()
    {
        return mCategory;
    }

    protected static String coverageString(final int coverage) { return format("Percentage%dX", coverage); }

    @Override
    public List<String> displayFieldNames() { return mFieldNames; }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            BamMetricSummary metrics = loadBamMetricsSummary(sampleId, fileSources.TumorBamMetrics);
            comparableItems.add(new BamMetricsData(mCategory, metrics, mComparisonPercentages, mFields));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load tumor BAM metrics data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }
}
