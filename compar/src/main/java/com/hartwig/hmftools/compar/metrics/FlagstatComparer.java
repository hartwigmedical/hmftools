package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.determineFlagStatsFilePath;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class FlagstatComparer extends ItemComparer
{
    public final CategoryType mCategory;

    protected enum Fields
    {
        MappedProportion;
    }

    public FlagstatComparer(final CategoryType category, final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);
        mCategory = category;

        mFields.add(new FieldInfo(
                Fields.MappedProportion.toString(),
                getOrMakeFieldCheck(
                        fieldCheckMap, Fields.MappedProportion.toString(), MAPPED_PROPORTION_ABS_THRESHOLD, MAPPED_PROPORTION_PCT_THRESHOLD),
                "%.2f"));
    }

    @Override
    public CategoryType category()
    {
        return mCategory;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            String sourceFile = isTumor() ? fileSources.TumorFlagstat : fileSources.GermlineFlagstat;
            String sourceSampleId = isTumor() ? sampleId : germlineSampleId;
            BamFlagStats flagstat = BamFlagStats.read(determineFlagStatsFilePath(sourceSampleId, sourceFile));
            comparableItems.add(new FlagstatData(mCategory, flagstat, mFields));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load {} flagstat data: {}", sampleId, isTumor() ? "tumor" : "germline", e.toString());
            return null;
        }
        return comparableItems;
    }

    private boolean isTumor() { return mCategory == CategoryType.TUMOR_FLAGSTAT; }
}
