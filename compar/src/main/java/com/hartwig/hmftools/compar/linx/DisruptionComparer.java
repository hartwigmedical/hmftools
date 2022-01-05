package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DisruptionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DisruptionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        List<LinxBreakend> breakends = dbAccess.readBreakends(sampleId);
        Map<String,DisruptionData> disruptions = Maps.newHashMap();

        breakends.forEach(x -> processBreakend(x, disruptions));

        return disruptions.values().stream().collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(fileSources.Linx, sampleId));
            Map<String,DisruptionData> disruptions = Maps.newHashMap();

            breakends.forEach(x -> processBreakend(x, disruptions));

            return disruptions.values().stream().collect(Collectors.toList());
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx breakend disruption data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }

    private void processBreakend(final LinxBreakend breakend, final Map<String,DisruptionData> disruptions)
    {
        if(!mConfig.DriverGenes.isEmpty() && !mConfig.DriverGenes.contains(breakend.gene()))
            return;

        String mappedName = mConfig.getGeneMappedName(breakend.gene());

        if(mappedName == null)
            return; // ignore genes that have been dropped

        DisruptionData disruptionData = disruptions.get(mappedName);

        if(disruptionData == null)
        {
            disruptionData = new DisruptionData(mappedName);
            disruptions.put(mappedName, disruptionData);
        }

        disruptionData.breakends().add(breakend);
    }

}
