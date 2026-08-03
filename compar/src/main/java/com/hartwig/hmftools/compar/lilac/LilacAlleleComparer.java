package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_ALLELE;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class LilacAlleleComparer extends ItemComparer
{
    protected enum Fields
    {
        RefTotal,
        TumorTotal,
        TumorCopyNumber,
        SomaticMissense,
        SomaticNonsenseOrFrameshift,
        SomaticSplice,
        SomaticInframeIndel,
        SomaticSynonymous;
    }

    private static final double FRAG_DIFF_PERC = 0.01;
    private static final double FRAG_DIFF_ABS = 10;
    private static final double VARIANT_DIFF_PERC = 0.1;
    private static final double VARIANT_DIFF_ABS = 0.4;

    public LilacAlleleComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.RefTotal.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.RefTotal.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.TumorTotal.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TumorTotal.toString(), FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                null));

        mFields.add(new FieldInfo(
                Fields.TumorCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TumorCopyNumber.toString(), 0.5, 0.15),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.SomaticMissense.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SomaticMissense.toString(), VARIANT_DIFF_ABS, VARIANT_DIFF_PERC),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.SomaticNonsenseOrFrameshift.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SomaticNonsenseOrFrameshift.toString(), VARIANT_DIFF_ABS, VARIANT_DIFF_PERC),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.SomaticSplice.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SomaticSplice.toString(), VARIANT_DIFF_ABS, VARIANT_DIFF_PERC),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.SomaticInframeIndel.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SomaticInframeIndel.toString(), VARIANT_DIFF_ABS, VARIANT_DIFF_PERC),
                "%.1f"));

        mFields.add(new FieldInfo(
                Fields.SomaticSynonymous.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SomaticSynonymous.toString(), VARIANT_DIFF_ABS, VARIANT_DIFF_PERC),
                "%.1f"));
    }

    @Override
    public CategoryType category() { return LILAC_ALLELE; }

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
            Map<String, Integer> alleleToSeen = new HashMap<>();
            for(LilacAllele allele : LilacAllele.read(LilacAllele.generateFilename(fileSources.Lilac, sampleId)))
            {
                String alleleString = String.format("%s:%s", allele.genes(), allele.allele());
                int index = alleleToSeen.getOrDefault(alleleString, 1);
                comparableItems.add(new LilacAlleleData(allele, index, mFields));

                alleleToSeen.put(alleleString, index + 1);
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Lilac allele data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
