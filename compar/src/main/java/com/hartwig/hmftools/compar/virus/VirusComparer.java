package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class VirusComparer extends ItemComparer
{
    protected enum Fields
    {
        Reported,
        Integrations,
        MeanCoverage,
        DriverLikelihood;
    }

    public VirusComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Reported.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Reported.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Integrations.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Integrations.toString(), null, 0.2),
                null));

        mFields.add(new FieldInfo(
                Fields.MeanCoverage.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MeanCoverage.toString(), null, 0.15),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.DriverLikelihood.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DriverLikelihood.toString()),null));
    }

    @Override
    public CategoryType category()
    {
        return VIRUS;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            AnnotatedVirusFile.read(AnnotatedVirusFile.generateFileName(fileSources.Virus, sampleId))
                    .forEach(v -> comparableItems.add(new VirusData(v, mFields)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Virus interpreter data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

}
