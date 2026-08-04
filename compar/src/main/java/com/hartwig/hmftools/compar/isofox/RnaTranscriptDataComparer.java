package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.TranscriptExpressionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaTranscriptDataComparer extends ItemComparer
{
    protected enum Fields
    {
        AdjTPM;
    }

    public RnaTranscriptDataComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.AdjTPM.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.AdjTPM.toString(), null, 0.05),
                "%.4e3"));
    }
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_TRANSCRIPT_DATA;
    }

    @Override
    public boolean hasReportable()
    {
        return false;
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
            String filename = determineFileName(sampleId, fileSources);
            TranscriptExpressionFile.read(filename).stream().map(x -> new RnaTranscriptData(x, mFields)).forEach(comparableItems::add);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Transcript data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final PipelineSourcePaths fileSources)
    {
        String filename = TranscriptExpressionFile.generateFilename(fileSources.Isofox, sampleId);
        String oldFilename = filename.replace(".tsv", ".csv");

        if(!Files.exists(Paths.get(filename)) && Files.exists(Paths.get(oldFilename)))
        {
            return oldFilename;
        }
        else
        {
            return filename;
        }
    }
}
