package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeTargetedRegions;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.mergeOverlapAndAdjacentRegions;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Holds the panel output data, including probes, target regions, and rejected features.
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
    public Stream<ChrBaseRegion> coveredRegions()
    {
        return mData.probes().stream().flatMap(probe -> probe.definition().regions().stream());
    }

    public void addResult(final ProbeGenerationResult result)
    {
        result.probes().forEach(probe ->
        {
            if(!probe.accepted())
            {
                throw new IllegalArgumentException("Should only add accepted probes to the panel");
            }
        });
        LOGGER.debug("Adding to panel: probes={} candidateTargetRegions={} rejectedFeatures={}",
                result.probes().size(), result.candidateTargetRegions().size(), result.rejectedFeatures().size());
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
        // Merge adjacent/overlapping target regions which have the same metadata.
        // If multiple target regions with different metadata overlap, there will be overlapping output regions.
        return mData.probes().stream()
                .collect(Collectors.groupingBy(Probe::metadata)).entrySet().stream()
                .flatMap(entry ->
                        mergeOverlapAndAdjacentRegions(entry.getValue().stream()
                                .flatMap(probe -> probeTargetedRegions(probe.definition(), probe.targetedRange()).stream()))
                                .stream()
                                .map(region -> new TargetRegion(region, entry.getKey())))
                .toList();
    }

    public List<RejectedFeature> rejectedFeatures()
    {
        return mData.rejectedFeatures();
    }
}
