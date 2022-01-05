package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

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
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
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
            CMP_LOGGER.info("sample({}) failed to load Linx fusion data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }

    private List<ComparableItem> processFusions(final List<LinxFusion> fusions)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        for(LinxFusion fusion : fusions)
        {
            if(fusion.reportedType().equals(KnownFusionType.NONE.toString()))
                continue;

            String mappedFusionName = getGeneMappedName(fusion.name());

            if(mappedFusionName == null)
                continue;

            comparableItems.add(new FusionData(fusion, mappedFusionName));
        }

        return comparableItems;
    }

    private String getGeneMappedName(final String fusionName)
    {
        String[] genes = fusionName.split("_");
        if(genes.length != 2)
            return fusionName;

        String upGeneName = mConfig.getGeneMappedName(genes[0]);
        String downGeneName = mConfig.getGeneMappedName(genes[1]);

        if(upGeneName == null || downGeneName == null)
            return null;

        return upGeneName  + "_" + downGeneName;
    }
}
