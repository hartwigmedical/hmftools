package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Holds the panel output data, including probes, target regions, and rejected regions.
// It's a mutable data structure because it's also used during probe generation to check which regions are already covered.
public class PanelData
{
    private ProbeGenerationResult mData;

    private static final Logger LOGGER = LogManager.getLogger(PanelData.class);

    public PanelData()
    {
        mData = new ProbeGenerationResult();
    }

    // Checks if a target region is fully covered by probes in the panel.
    public boolean isCovered(final ChrBaseRegion target)
    {
        return probeRegions().anyMatch(probe -> probe.containsRegion(target));
    }

    private Stream<ChrBaseRegion> probeRegions()
    {
        return mData.probes().stream().map(probe -> probe.candidate().probeRegion());
    }

    public void addResult(final ProbeGenerationResult result)
    {
        LOGGER.debug("Adding to panel: probes={} targetRegions={} rejectedRegions={}",
                result.probes().size(), result.targetRegions().size(), result.rejectedRegions().size());
        mData = mData.add(result);
    }

    public List<EvaluatedProbe> probes()
    {
        return mData.probes();
    }

    public List<TargetRegion> targetRegions()
    {
        return mData.targetRegions();
    }

    public List<RejectedRegion> rejectedRegions()
    {
        return mData.rejectedRegions();
    }
}
