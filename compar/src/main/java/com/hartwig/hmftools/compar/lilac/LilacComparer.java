package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_MISSENSE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_NFS;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SPLICE;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.common.Category.LILAC;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_ALLELES;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_VARIANTS;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SampleFileSources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class LilacComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    private static final double FRAG_DIFF_PERC = 0.01;
    private static final double FRAG_DIFF_ABS = 10;
    private static final double VARIANT_DIFF_PERC = 0.1;
    private static final double VARIANT_DIFF_ABS = 0.4;

    public LilacComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return LILAC; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(LilacAllele.FLD_REF_TOTAL, FRAG_DIFF_ABS, FRAG_DIFF_PERC);
        thresholds.addFieldThreshold(LilacAllele.FLD_TUMOR_TOTAL, FRAG_DIFF_ABS, FRAG_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_TOTAL_FRAGS, FRAG_DIFF_ABS, FRAG_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_FIT_FRAGS, FRAG_DIFF_ABS, FRAG_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_DISC_ALIGN_FRAGS, FRAG_DIFF_ABS, FRAG_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_DISC_INDELS, FRAG_DIFF_ABS, FRAG_DIFF_PERC);

        thresholds.addFieldThreshold(FLD_MISSENSE, VARIANT_DIFF_ABS, VARIANT_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_NFS, VARIANT_DIFF_ABS, VARIANT_DIFF_PERC);
        thresholds.addFieldThreshold(FLD_SPLICE, VARIANT_DIFF_ABS, VARIANT_DIFF_PERC);
        thresholds.addFieldThreshold(LilacAllele.FLD_INDEL, VARIANT_DIFF_ABS, VARIANT_DIFF_PERC);
        thresholds.addFieldThreshold(LilacAllele.FLD_SYNON, VARIANT_DIFF_ABS, VARIANT_DIFF_PERC);
        thresholds.addFieldThreshold(LilacAllele.FLD_TUMOR_CN, 0.5, 0.15);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_QC_STATUS, FLD_ALLELES, FLD_VARIANTS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // Not currently supported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final SampleFileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            LilacQcData qcData = LilacQcData.read(LilacQcData.generateFilename(fileSources.lilac(), sampleId));
            List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(fileSources.lilac(), sampleId));

            comparableItems.add(new LilacData(qcData, alleles));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Lilac data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
