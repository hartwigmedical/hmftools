package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.isPositionValid;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSvProbe;

import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// TODO: support single breakends?

// Probes covering a list of arbitrary structural variants provided by the user.
public class CustomSvs
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CUSTOM_SV;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CUSTOM_SV_QUALITY_MIN, CUSTOM_SV_GC_TARGET, CUSTOM_SV_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(CustomSvs.class);

    public static void generateProbes(final String customSvFile, final Map<String, Integer> chromosomeLengths,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating custom structural variant probes");

        List<CustomSv> customSvs = CustomSv.readFromFile(customSvFile);

        checkSvPositions(customSvs, chromosomeLengths);
        checkNoDuplicates(customSvs);

        ProbeGenerationResult result = customSvs.stream()
                .map(customSv -> generateProbes(customSv, probeGenerator, panelData))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
        // TODO: should try to check overlaps better?
        // Generate all the probes at once and then add to the result because it's too hard to check overlap in advance, and it's unlikely
        // the user accidentally specified overlapping SVs.
        panelData.addResult(result);

        LOGGER.info("Done generating custom structural variant probes");
    }

    private static void checkSvPositions(final List<CustomSv> customSvs, final Map<String, Integer> chromosomeLengths)
    {
        List<CustomSv> invalid = customSvs.stream()
                .filter(sv ->
                        !(isPositionValid(sv.startPosition(), chromosomeLengths) && isPositionValid(sv.endPosition(), chromosomeLengths)))
                .toList();
        if(!invalid.isEmpty())
        {
            invalid.forEach(sv -> LOGGER.error("Invalid custom structural variant positions: {}", sv));
            throw new UserInputError("Invalid custom structural variant positions");
        }
    }

    private static void checkNoDuplicates(final List<CustomSv> customSvs)
    {
        LOGGER.debug("Checking custom structural variants for duplicates");
        List<CustomSv> duplicated = customSvs.stream()
                .filter(sv ->
                        customSvs.stream().anyMatch(sv2 -> sv != sv2
                                && sv.startPosition() == sv2.startPosition() && sv.startOrientation() == sv2.startOrientation()
                                && sv.endPosition() == sv2.endPosition() && sv.endOrientation() == sv2.endOrientation())
                ).toList();
        if(!duplicated.isEmpty())
        {
            duplicated.forEach(sv -> LOGGER.error("Duplicate custom structural variant: {}", sv));
            throw new UserInputError("Duplicate custom structural variants");
        }
    }

    private static ProbeGenerationResult generateProbes(final CustomSv sv, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        LOGGER.debug("Generating probes for {}", sv);
        TargetMetadata metadata = new TargetMetadata(TARGET_TYPE, sv.extraInfo());
        SequenceDefinition definition = buildSvProbe(
                sv.startPosition().Chromosome, sv.startPosition().Position, sv.startOrientation(),
                sv.endPosition().Chromosome, sv.endPosition().Position, sv.endOrientation(),
                sv.insertSequence(),
                PROBE_LENGTH);
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        ProbeGenerationResult result = probeGenerator.probe(definition, targetedRange, metadata, PROBE_CRITERIA, coverage);
        return result;
    }
}

