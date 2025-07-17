package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CUSTOM_REGION_QUALITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TOLERANCE;

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
    private static final TargetRegionType TARGET_REGION_TYPE = TargetRegionType.CUSTOM;

    private static final ProbeSelectCriteria PROBE_SELECT_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvalCriteria(CUSTOM_REGION_QUALITY_MIN, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE),
            new ProbeSelectStrategy.MaxQuality());

    private static final String FLD_EXTRA_INFO = "ExtraInfo";

    private static final Logger LOGGER = LogManager.getLogger(CustomRegions.class);

    public static ProbeGenerationResult generateProbes(final String customRegionFile, final Map<String, Integer> chromosomeLengths,
            final ProbeGenerator probeGenerator)
    {
        LOGGER.info("Generating custom region probes");

        List<CustomRegion> customRegions = loadCustomRegionsFile(customRegionFile);

        checkRegionBounds(customRegions, chromosomeLengths);
        checkNoOverlaps(customRegions);

        // TODO: handle overlaps with existing probes

        ProbeGenerationResult result = customRegions.stream()
                .map(region -> generateProbes(region, probeGenerator))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        LOGGER.info("Done generating custom region probes");
        return result;
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
            int chromosomeIdx = reader.getColumnIndex(FLD_CHROMOSOME);
            int posStartIdx = reader.getColumnIndex(FLD_POSITION_START);
            int posEndIdx = reader.getColumnIndex(FLD_POSITION_END);
            int extraInfoIdx = reader.getColumnIndex(FLD_EXTRA_INFO);

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

    private static ProbeGenerationResult generateProbes(final CustomRegion region, final ProbeGenerator probeGenerator)
    {
        LOGGER.debug("Generating probes for {}", region);
        TargetMetadata metadata = new TargetMetadata(TARGET_REGION_TYPE, region.extraInfo());
        TargetRegion target = new TargetRegion(region.region(), metadata);
        ProbeGenerationResult result = probeGenerator.coverRegion(target, PROBE_SELECT_CRITERIA);
        return result;
    }
}
