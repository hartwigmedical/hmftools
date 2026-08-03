package com.hartwig.hmftools.compar.peach;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.PEACH;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class PeachComparer extends ItemComparer
{
    protected enum Fields
    {
        AlleleCount,
        Function,
        Drugs,
        PrescriptionUrls;
    }

    public PeachComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.AlleleCount.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.AlleleCount.toString(), null, null),
                null));

        mFields.add(new FieldInfo(
                Fields.Function.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Function.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Drugs.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Drugs.toString()), null));

        mFields.add(new FieldInfo(
                Fields.PrescriptionUrls.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.PrescriptionUrls.toString()), null));
    }

    @Override
    public CategoryType category()
    {
        return PEACH;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        String fileName = determineFileName(sampleId, germlineSampleId, fileSources);
        try
        {
            PeachGenotypeFile.read(fileName).forEach(g -> comparableItems.add(new PeachData(g, mFields)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Peach data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final String currentFileName = PeachGenotypeFile.generateFileName(fileSources.Peach, germlineSampleId);
        final String oldFileName = PeachGenotypeFile.generateOldPythonFileName(fileSources.Peach, sampleId);
        if(!fileExists(currentFileName) && fileExists(oldFileName))
        {
            return oldFileName;
        }
        else
        {
            return currentFileName;
        }
    }
}
