package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CN_SEGS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CONTAMINATION;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_FIT_METHOD;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GENDER;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GERM_ABS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_INDELS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PLOIDY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PURITY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_QC_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_SV_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_UNS_CN_SEGS;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class PurityComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public PurityComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return PURITY; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_PURITY, 0.01, 0);
        thresholds.addFieldThreshold(FLD_PLOIDY, 0.1, 0);
        thresholds.addFieldThreshold(FLD_CONTAMINATION, 0.005, 0);
        thresholds.addFieldThreshold(FLD_TMB, 0, 0.01);
        thresholds.addFieldThreshold(FLD_MS_INDELS, 0, 0.01);
        thresholds.addFieldThreshold(FLD_TML, 0, 0.01);
        thresholds.addFieldThreshold(FLD_CN_SEGS, 0, 0.1);
        thresholds.addFieldThreshold(FLD_UNS_CN_SEGS, 0, 0.1);
        thresholds.addFieldThreshold(FLD_SV_TMB, 0, 0.03);
    }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_PURITY, FLD_PLOIDY, FLD_CONTAMINATION, FLD_TMB, FLD_TML, FLD_MS_INDELS, FLD_SV_TMB, FLD_CN_SEGS ,FLD_UNS_CN_SEGS,
                FLD_QC_STATUS, FLD_GENDER, FLD_GERM_ABS, FLD_FIT_METHOD, FLD_MS_STATUS, FLD_TMB_STATUS, FLD_TML_STATUS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final PurityContext purityContext = dbAccess.readPurityContext(sampleId);

        List<ComparableItem> items = Lists.newArrayList();
        items.add(new PurityData(purityContext));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            PurityContext purityContext = PurityContextFile.read(fileSources.Purple, sampleId);
            comparableItems.add(new PurityData(purityContext));

        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Purple purity data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }
}
