package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_QC;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class LilacQcComparer extends ItemComparer
{
    protected enum Fields
    {
        Status,
        TotalFragments,
        FittedFragments,
        DiscardedAlignmentFragments,
        DiscardedIndels,
        HlaYAllele,
        Alleles;
    }

    private static final double FRAG_DIFF_PERC = 0.01;
    private static final double FRAG_DIFF_ABS = 10;

    public LilacQcComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Status.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Status.toString()), null));

        mFields.add(new FieldInfo(
                Fields.TotalFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TotalFragments.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.FittedFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FittedFragments.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.DiscardedAlignmentFragments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DiscardedAlignmentFragments.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.DiscardedIndels.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DiscardedIndels.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.HlaYAllele.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.HlaYAllele.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Alleles.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Alleles.toString()), null));
    }

    @Override
    public CategoryType category() { return LILAC_QC; }

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
            // add an item for each gene class
            List<LilacQcData> qcDataList = LilacQcData.read(LilacQcData.generateFilename(fileSources.Lilac, sampleId));
            List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(fileSources.Lilac, sampleId));

            for(LilacQcData qcData : qcDataList)
            {
                List<LilacAllele> geneAlleles = alleles.stream()
                        .filter(x -> x.genes().equals(qcData.genes()))
                        .collect(Collectors.toList());
                comparableItems.add(new LilacQcComparData(qcData, geneAlleles, mFields));
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Lilac data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
