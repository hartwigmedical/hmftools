package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_INDEL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_MISSENSE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_NFS;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_REF_TOTAL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SPLICE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SYNON;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_CN;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_TOTAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_ALLELE;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class LilacAlleleComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    private static final double FRAG_DIFF_PERC = 0.01;
    private static final double FRAG_DIFF_ABS = 10;
    private static final double VARIANT_DIFF_PERC = 0.1;
    private static final double VARIANT_DIFF_ABS = 0.4;

    public LilacAlleleComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return LILAC_ALLELE; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        boolean includeDetailedFields = matchLevel.equals(MatchLevel.DETAILED);
        return List.of(
                new DoubleField(FLD_MISSENSE, i -> ((LilacAlleleData) i).Allele.somaticMissense(), true,
                        VARIANT_DIFF_ABS, VARIANT_DIFF_PERC, "%.2f"),
                new DoubleField(FLD_NFS, i -> ((LilacAlleleData) i).Allele.somaticNonsenseOrFrameshift(), true,
                        VARIANT_DIFF_ABS, VARIANT_DIFF_PERC, "%.2f"),
                new DoubleField(FLD_SPLICE, i -> ((LilacAlleleData) i).Allele.somaticSplice(), true,
                        VARIANT_DIFF_ABS, VARIANT_DIFF_PERC, "%.2f"),
                new DoubleField(FLD_INDEL, i -> ((LilacAlleleData) i).Allele.somaticInframeIndel(), true,
                        VARIANT_DIFF_ABS, VARIANT_DIFF_PERC, "%.2f"),
                new DoubleField(FLD_TUMOR_CN, i -> ((LilacAlleleData) i).Allele.tumorCopyNumber(), true,
                        0.5, 0.15, "%.2f"),
                new IntField(FLD_REF_TOTAL, i -> ((LilacAlleleData) i).Allele.refFragments(), includeDetailedFields,
                        FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new IntField(FLD_TUMOR_TOTAL, i -> ((LilacAlleleData) i).Allele.tumorFragments(), includeDetailedFields,
                        FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new DoubleField(FLD_SYNON, i -> ((LilacAlleleData) i).Allele.somaticSynonymous(), includeDetailedFields,
                        VARIANT_DIFF_ABS, VARIANT_DIFF_PERC, "%.2f")
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        // Not currently supported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            Map<String, Integer> alleleToSeen = new HashMap<>();
            for(LilacAllele allele : LilacAllele.read(LilacAllele.generateFilename(fileSources.Lilac, sampleId)))
            {
                String alleleString = String.format("%s:%s", allele.genes(), allele.allele());
                int index = alleleToSeen.getOrDefault(alleleString, 1);
                comparableItems.add(new LilacAlleleData(allele, index));

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
