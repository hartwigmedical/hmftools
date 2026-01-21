package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_GERMLINE_CN;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_GERMLINE_STATUS;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_TUMOR_CN;
import static com.hartwig.hmftools.compar.purple.GermlineAmpDelData.FLD_TUMOR_STATUS;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineAmpDelComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineAmpDelComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return GERMLINE_AMP_DEL; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_GERMLINE_CN, 0.2, 0.1);
        thresholds.addFieldThreshold(FLD_TUMOR_CN, 0.2, 0.1);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_GERMLINE_STATUS, FLD_TUMOR_STATUS, FLD_GERMLINE_CN, FLD_TUMOR_CN);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<GermlineAmpDel> germlineAmpDels = dbAccess.readGermlineCopyNumbers(sampleId);
        List<ComparableItem> items = Lists.newArrayList();
        germlineAmpDels.forEach(x -> items.add(createGermlineAmpDelData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            String germlineAmpDelFile = GermlineAmpDel.generateFilename(fileSources.Purple, sampleId);

            if(!Files.exists(Paths.get(germlineAmpDelFile)))
            {
                // try pre v3.0 germline deletions
                germlineAmpDelFile = germlineAmpDelFile.replaceAll("purple.germline_amp_del.tsv", "purple.germline.deletion.tsv");
            }

            List<GermlineAmpDel> germlineAmpDels = GermlineAmpDel.read(germlineAmpDelFile);
            germlineAmpDels.forEach(x -> comparableItems.add(createGermlineAmpDelData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read germline amp-del data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private GermlineAmpDelData createGermlineAmpDelData(final GermlineAmpDel deletion)
    {
        String comparisonChromosome = determineComparisonChromosome(deletion.Chromosome, mConfig.RequiresLiftover);
        return new GermlineAmpDelData(deletion, comparisonChromosome);
    }
}
