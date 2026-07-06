package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_HLA_Y;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_ALLELES;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.compar.common.field.StringListField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class LilacComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    private static final double FRAG_DIFF_PERC = 0.01;
    private static final double FRAG_DIFF_ABS = 10;

    public LilacComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return LILAC; }

    @Override
    public List<Field> fields()
    {
        return List.of(
                new StringField(FLD_QC_STATUS, i -> ((LilacData) i).QcData.status(), true),
                new IntField(FLD_TOTAL_FRAGS, i -> ((LilacData) i).QcData.totalFragments(), true, FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new IntField(FLD_FIT_FRAGS, i -> ((LilacData) i).QcData.fittedFragments(), true, FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new IntField(FLD_DISC_ALIGN_FRAGS, i -> ((LilacData) i).QcData.discardedAlignmentFragments(), true, FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new IntField(FLD_DISC_INDELS, i -> ((LilacData) i).QcData.discardedIndels(), true, FRAG_DIFF_ABS, FRAG_DIFF_PERC),
                new StringField(FLD_HLA_Y, i -> ((LilacData) i).QcData.hlaYAllele(), true),
                new StringListField(FLD_ALLELES, i -> ((LilacData) i).Alleles.stream()
                        .map(LilacAllele::allele).sorted().toList(), true)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(FLD_QC_STATUS, FLD_ALLELES);
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
            // add an item for each gene class
            List<LilacQcData> qcDataList = LilacQcData.read(LilacQcData.generateFilename(fileSources.Lilac, sampleId));
            List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(fileSources.Lilac, sampleId));

            for(LilacQcData qcData : qcDataList)
            {
                List<LilacAllele> geneAlleles = alleles.stream().filter(x -> x.genes().equals(qcData.genes())).collect(Collectors.toList());
                comparableItems.add(new LilacData(qcData, geneAlleles));
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
