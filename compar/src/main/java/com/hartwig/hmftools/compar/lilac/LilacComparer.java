package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.Category.LILAC;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_ALLELES;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_STATUS;
import static com.hartwig.hmftools.compar.lilac.LilacData.FLD_VARIANTS;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class LilacComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public LilacComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return LILAC; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_STATUS, FLD_ALLELES, FLD_VARIANTS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            LilacSummaryData lilacData = LilacSummaryData.read(fileSources.Lilac, sampleId);

            comparableItems.add(new LilacData(lilacData));
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Lilac data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }
}
