package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_REGION_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_REGION_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_REGION_QUALITY_MIN_DEFAULT;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.isRegionValid;
import static com.hartwig.hmftools.panelbuilder.Utils.findDuplicates;

import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering a list of arbitrary regions provided by the user.
public class CustomRegions
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CUSTOM_REGION;

    private static final ProbeSelector.Strategy PROBE_SELECT = new ProbeSelector.Strategy.MaxQuality();

    private static final Logger LOGGER = LogManager.getLogger(CustomRegions.class);

    public static void generateProbes(final List<String> customRegionFiles, final Map<String, Integer> chromosomeLengths,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating custom region probes");

        List<CustomRegion> customRegions = customRegionFiles.stream()
                .flatMap(customRegionFile -> CustomRegion.readFromFile(customRegionFile).stream())
                .toList();

        checkRegionBounds(customRegions, chromosomeLengths);
        checkNoOverlaps(customRegions);

        ProbeGenerationResult result = generateProbes(customRegions, probeGenerator, panelData);
        // Mostly safe to generate the probes at once and then add to the result afterward,
        // since we already checked that the custom regions don't overlap and assume any overlap from the probe generation is minimal.
        // Potentially, there could be a small overlap if two custom regions were right next to each other.
        panelData.addResult(result);

        LOGGER.info("Done generating custom region probes");
    }

    private static void checkRegionBounds(final List<CustomRegion> customRegions, final Map<String, Integer> chromosomeLengths)
    {
        List<CustomRegion> invalid = customRegions.stream()
                .filter(region -> !isRegionValid(region.region(), chromosomeLengths))
                .toList();
        if(!invalid.isEmpty())
        {
            invalid.forEach(region -> LOGGER.error("Invalid custom region bounds: {}", region));
            throw new UserInputError("Invalid custom region bounds");
        }
    }

    private static void checkNoOverlaps(final List<CustomRegion> customRegions)
    {
        LOGGER.debug("Checking custom regions for overlap");
        List<CustomRegion> overlapped = findDuplicates(customRegions, (region1, region2) -> region1.region().overlaps(region2.region()));
        if(!overlapped.isEmpty())
        {
            overlapped.forEach(region -> LOGGER.error("Overlapping custom region: {}", region));
            throw new UserInputError("Custom regions overlap");
        }
    }

    private static ProbeGenerationResult generateProbes(final List<CustomRegion> customRegions, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        Stream<ProbeGenerationSpec> probeGenerationSpecs = customRegions.stream().map(CustomRegions::createProbeGenerationSpec);
        return probeGenerator.generateBatch(probeGenerationSpecs, coverage);
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final CustomRegion region)
    {
        LOGGER.debug("Generating probes for {}", region);
        TargetMetadata metadata = new TargetMetadata(TARGET_TYPE, region.extraInfo());
        ProbeEvaluator.Criteria evalCriteria = new ProbeEvaluator.Criteria(
                region.qualityScoreMin() == null ? CUSTOM_REGION_QUALITY_MIN_DEFAULT : region.qualityScoreMin(),
                CUSTOM_REGION_GC_TARGET, CUSTOM_REGION_GC_TOLERANCE);
        return new ProbeGenerationSpec.CoverRegion(region.region(), metadata, evalCriteria, PROBE_SELECT);
    }
}
