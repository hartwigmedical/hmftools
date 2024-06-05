package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.Category.FUSION;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_LINKS;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_TERM;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_KEPT;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_LOST;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_DOWN;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_UP;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_PHASED;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_REPORTED_TYPE;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class FusionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public FusionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return FUSION; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_REPORTED_TYPE, FLD_PHASED, FLD_LIKELIHOOD, FLD_EXON_UP, FLD_EXON_DOWN, FLD_CHAIN_LINKS, FLD_CHAIN_TERM,
                FLD_DOMAINS_KEPT, FLD_DOMAINS_LOST);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        List<LinxFusion> fusions = dbAccess.readFusions(sampleId);
        return processFusions(fusions);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(fileSources.Linx, sampleId));

            return processFusions(fusions);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx fusion data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> processFusions(final List<LinxFusion> fusions)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        for(LinxFusion fusion : fusions)
        {
            if(fusion.reportedType().equals(KnownFusionType.NONE.toString()))
                continue;

            comparableItems.add(new FusionData(fusion, fusion.name()));
        }

        return comparableItems;
    }
}
