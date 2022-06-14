package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.GERMLINE_DELETION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_GERMLINE_CN;
import static com.hartwig.hmftools.compar.purple.GermlineDeletionData.FLD_TUMOR_CN;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineDeletionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineDeletionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return GERMLINE_DELETION; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_GERMLINE_CN, 0.2, 0.1);
        thresholds.addFieldThreshold(FLD_TUMOR_CN, 0.2, 0.1);
    }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<GermlineDeletion> germlineDeletions = dbAccess.readGermlineDeletions(sampleId);
        List<ComparableItem> items = Lists.newArrayList();
        germlineDeletions.forEach(x -> items.add(new GermlineDeletionData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<GermlineDeletion> germlineDeletions = GermlineDeletion.read(GermlineDeletion.generateFilename(fileSources.Purple, sampleId));
            germlineDeletions.forEach(x -> comparableItems.add(new GermlineDeletionData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to read germline deletion data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }
}
