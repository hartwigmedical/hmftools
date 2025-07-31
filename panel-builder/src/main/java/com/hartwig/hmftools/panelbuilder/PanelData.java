package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.Utils.isCoveredBy;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Holds the panel output data, including probes, target regions, and rejected regions.
// It's a mutable data structure because it's also used during probe generation to check which regions are already covered.
public class PanelData implements PanelCoverage
{
    private ProbeGenerationResult mData;

    private static final Logger LOGGER = LogManager.getLogger(PanelData.class);

    public PanelData()
    {
        mData = new ProbeGenerationResult();
    }

    @Override
    public boolean isCovered(final ChrBaseRegion region)
    {
        return isCoveredBy(region, coveredRegions());
    }

    @Override
    public Stream<ChrBaseRegion> coveredRegions()
    {
        return mData.probes().stream().map(Probe::region);
    }

    public void addResult(final ProbeGenerationResult result)
    {
        LOGGER.debug("Adding to panel: probes={} candidateTargetRegions={} coveredTargetRegions={} rejectedRegions={}",
                result.probes().size(), result.candidateTargetRegions().size(), result.coveredTargetRegions().size(),
                result.rejectedRegions().size());
        mData = mData.add(result);
    }

    public List<Probe> probes()
    {
        return mData.probes();
    }

    // All the target regions, regardless of how they were covered by probes.
    // This is mostly useful for informational and debugging purposes, since the regions may or may not be meaningful depending on the probe
    // source type.
    public List<TargetRegion> candidateTargetRegions()
    {
        return mData.candidateTargetRegions();
    }

    // The target regions which the probes aim to hit. This is the intersection of the probe and its target region.
    public List<TargetRegion> coveredTargetRegions()
    {
        return mData.coveredTargetRegions();
    }

    public List<RejectedRegion> rejectedRegions()
    {
        return mData.rejectedRegions();
    }
}
