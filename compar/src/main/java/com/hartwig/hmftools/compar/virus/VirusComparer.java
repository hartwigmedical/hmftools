package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.viridian.VirusReportFile;
import com.hartwig.hmftools.common.virusinterpreter.AnnotatedVirus;
import com.hartwig.hmftools.common.virusinterpreter.AnnotatedVirusFile;
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
        Name,
        Reported,
        Interpretation,
        Integrations,
        MeanCoverage,
        PercentageCovered,
        CopiesPerTumorCell,
        ClonalityStatus;
    }

    public VirusComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(Fields.Name.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Name.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Reported.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Reported.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Interpretation.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Interpretation.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Integrations.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Integrations.toString()), null));

        mFields.add(new FieldInfo(
                Fields.MeanCoverage.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MeanCoverage.toString(), null, 0.1),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.PercentageCovered.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.PercentageCovered.toString(), null, 0.1),
                "%.3f"));

        mFields.add(new FieldInfo(
                Fields.CopiesPerTumorCell.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.CopiesPerTumorCell.toString(), null, 0.1),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.ClonalityStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.ClonalityStatus.toString()), null));
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
        return fileSources.Viridian != null
                ? loadViridian(sampleId, fileSources.Viridian)
                : loadVirusInterpreter(sampleId, fileSources.VirusInterpreter);
    }

    private List<ComparableItem> loadViridian(final String sampleId, final String viridianDir)
    {
        try
        {
            return VirusReportFile.read(VirusReportFile.generateFilename(viridianDir, sampleId)).stream()
                    .map(x -> (ComparableItem) VirusData.fromViridian(x, mFields))
                    .collect(Collectors.toList());
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Viridian data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> loadVirusInterpreter(final String sampleId, final String virusInterpreterDir)
    {
        List<AnnotatedVirus> viruses;
        try
        {
            viruses = AnnotatedVirusFile.read(AnnotatedVirusFile.generateFileName(virusInterpreterDir, sampleId));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Virus interpreter data: {}", sampleId, e.toString());
            return null;
        }

        List<ComparableItem> items = viruses.stream()
                .map(x -> (ComparableItem) VirusData.fromVirusInterpreter(x, mFields))
                .collect(Collectors.toList());

        List<String> duplicateKeys = duplicateKeys(items);
        if(!duplicateKeys.isEmpty())
        {
            CMP_LOGGER.warn("sample({}) has more than one Virus interpreter row per virus group({})", sampleId, duplicateKeys);
            return null;
        }

        return items;
    }

    static List<String> duplicateKeys(final List<ComparableItem> items)
    {
        return items.stream()
                .collect(Collectors.groupingBy(ComparableItem::key, TreeMap::new, Collectors.counting()))
                .entrySet().stream()
                .filter(x -> x.getValue() > 1)
                .map(Map.Entry::getKey)
                .toList();
    }
}
