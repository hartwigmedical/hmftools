package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_REGION_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENERAL_GC_TOLERANCE;

import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering a list of arbitrary regions provided by the user.
public class CustomRegions
{
    private static final TargetMetadata.Type TARGET_REGION_TYPE = TargetMetadata.Type.CUSTOM;

    private static final ProbeSelectCriteria PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(CUSTOM_REGION_QUALITY_MIN, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE),
            new ProbeSelector.Strategy.MaxQuality());

    private static final String FLD_EXTRA_INFO = "ExtraInfo";

    private static final Logger LOGGER = LogManager.getLogger(CustomRegions.class);

    public static void generateProbes(final String customRegionFile, final Map<String, Integer> chromosomeLengths,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating custom region probes");

        List<CustomRegion> customRegions = loadCustomRegionsFile(customRegionFile);

        checkRegionBounds(customRegions, chromosomeLengths);
        checkNoOverlaps(customRegions);

        ProbeGenerationResult result = customRegions.stream()
                .map(region -> generateProbes(region, probeGenerator, panelData))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
        // Mostly safe to generate the probes at once and then add to the result afterward,
        // since we already checked that the custom region don't overlap, and assume any overlap from the probe generation is minimal.
        // Potentially there could be small overlap if two custom regions were right next to each other.
        panelData.addResult(result);

        LOGGER.info("Done generating custom region probes");
    }

    private record CustomRegion(
            ChrBaseRegion region,
            // Arbitrary descriptor for the user.
            String extraInfo
    )
    {
    }

    private static List<CustomRegion> loadCustomRegionsFile(final String filePath)
    {
        LOGGER.debug("Loading custom regions file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeIdx = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int posStartIdx = requireNonNull(reader.getColumnIndex(FLD_POSITION_START));
            int posEndIdx = requireNonNull(reader.getColumnIndex(FLD_POSITION_END));
            int extraInfoIdx = requireNonNull(reader.getColumnIndex(FLD_EXTRA_INFO));

            List<CustomRegion> regions = reader.stream().map(row ->
            {
                String chromosome = row.get(chromosomeIdx);
                int start = row.getInt(posStartIdx);
                int end = row.getInt(posEndIdx);
                String extraInfo = row.get(extraInfoIdx);
                ChrBaseRegion baseRegion = new ChrBaseRegion(chromosome, start, end);
                return new CustomRegion(baseRegion, extraInfo);
            }).toList();

            LOGGER.info("Loaded {} custom regions from {}", regions.size(), filePath);
            return regions;
        }
    }

    private static void checkRegionBounds(final List<CustomRegion> customRegions, final Map<String, Integer> chromosomeLengths)
    {
        Predicate<ChrBaseRegion> isRegionValid = region ->
        {
            if(!region.hasValidPositions())
            {
                return false;
            }
            Integer chromosomeLength = chromosomeLengths.get(region.chromosome());
            if(chromosomeLength == null)
            {
                return false;
            }
            return region.end() <= chromosomeLength;
        };

        List<CustomRegion> invalid = customRegions.stream()
                .filter(region -> !isRegionValid.test(region.region()))
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
        List<CustomRegion> invalid = customRegions.stream()
                .filter(region ->
                        customRegions.stream().anyMatch(region2 ->
                                region2 != region && region.region().overlaps(region2.region()))
                ).toList();
        if(!invalid.isEmpty())
        {
            invalid.forEach(region -> LOGGER.error("Overlapping custom region: {}", region));
            throw new UserInputError("Custom regions overlap");
        }
    }

    private static ProbeGenerationResult generateProbes(final CustomRegion region, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        LOGGER.debug("Generating probes for {}", region);
        TargetMetadata metadata = new TargetMetadata(TARGET_REGION_TYPE, region.extraInfo());
        ProbeGenerationResult result = probeGenerator.coverRegion(region.region(), metadata, PROBE_CRITERIA, coverage);
        return result;
    }
}
