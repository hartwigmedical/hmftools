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
        return dbAccess.readFusions(sampleId).stream()
                .filter(x -> !x.reportedType().equals(KnownFusionType.NONE.toString()))
                .map(x -> new FusionData(x)).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = com.google.common.collect.Lists.newArrayList();

        try
        {
            List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(fileSources.Linx, sampleId));

            fusions.stream()
                    .filter(x -> !x.reportedType().equals(KnownFusionType.NONE.toString()))
                    .forEach(x -> comparableItems.add(new FusionData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx fusion data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }
}
