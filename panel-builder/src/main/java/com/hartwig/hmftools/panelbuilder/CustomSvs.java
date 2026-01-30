package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SV_QUALITY_MIN_DEFAULT;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.isPositionValid;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSvProbe;
import static com.hartwig.hmftools.panelbuilder.Utils.findDuplicates;

import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// TODO? support single breakends?

// Probes covering a list of arbitrary structural variants provided by the user.
public class CustomSvs
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CUSTOM_SV;

    private static final Logger LOGGER = LogManager.getLogger(CustomSvs.class);

    public static void generateProbes(final String customSvFile, final Map<String, Integer> chromosomeLengths,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating custom structural variant probes");

        List<CustomSv> customSvs = CustomSv.readFromFile(customSvFile);

        checkSvPositions(customSvs, chromosomeLengths);
        checkNoDuplicates(customSvs);

        generateProbes(customSvs, probeGenerator, panelData);

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
        List<CustomSv> duplicated = findDuplicates(customSvs, (sv1, sv2) ->
                sv1.startPosition() == sv2.startPosition() && sv1.startOrientation() == sv2.startOrientation() &&
                        sv1.endPosition() == sv2.endPosition() && sv1.endOrientation() == sv2.endOrientation());
        if(!duplicated.isEmpty())
        {
            duplicated.forEach(customSv -> LOGGER.error("Duplicate custom structural variant: {}", customSv));
            throw new UserInputError("Duplicate custom structural variants");
        }
    }

    private static void generateProbes(final List<CustomSv> customSvs, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        Stream<ProbeGenerationSpec> probeGenerationSpecs = customSvs.stream().map(CustomSvs::createProbeGenerationSpec);
        probeGenerator.generateBatch(probeGenerationSpecs, panelData);
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final CustomSv customSv)
    {
        LOGGER.debug("Generating probes for {}", customSv);
        TargetMetadata metadata = new TargetMetadata(TARGET_TYPE, customSv.extraInfo());
        SequenceDefinition definition = buildSvProbe(
                customSv.startPosition().Chromosome, customSv.startPosition().Position, customSv.startOrientation(),
                customSv.endPosition().Chromosome, customSv.endPosition().Position, customSv.endOrientation(),
                customSv.insertSequence(),
                PROBE_LENGTH);
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        ProbeEvaluator.Criteria evalCriteria = new ProbeEvaluator.Criteria(
                customSv.qualityScoreMin() == null ? CUSTOM_SV_QUALITY_MIN_DEFAULT : customSv.qualityScoreMin(),
                CUSTOM_SV_GC_TARGET, CUSTOM_SV_GC_TOLERANCE);
        return new ProbeGenerationSpec.SingleProbe(definition, targetedRange, metadata, evalCriteria);
    }
}

