package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import org.jetbrains.annotations.NotNull;

public class RnaGeneDataComparer extends ItemComparer
{
    protected enum Fields
    {
        SplicedFragments,
        UnsplicedFragments,
        AdjTPM;
    }

    public RnaGeneDataComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.SplicedFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SplicedFragments.toString(), 10.0, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.UnsplicedFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.UnsplicedFragments.toString(), 10.0, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.AdjTPM.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.AdjTPM.toString(), null, 0.05),
                "%.4e3"));
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_GENE_DATA;
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
        String filename = GeneExpressionFile.generateFilename(fileSources.Isofox, sampleId);
        filename = checkOldIsofoxFilename(filename);

        List<GeneExpression> geneExpressions = GeneExpressionFile.read(filename);
        if(geneExpressions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Gene data", sampleId);
            return null;
        }

        return geneExpressions.stream().map(x -> new RnaGeneData(x, mFields)).collect(Collectors.toList());
    }

    protected static String checkOldIsofoxFilename(final String filename)
    {
        String oldFilename = filename.replace(TSV_EXTENSION, CSV_EXTENSION);

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
