package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SMALL_VARIANT_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SMALL_VARIANT_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CUSTOM_SMALL_VARIANT_QUALITY_MIN_DEFAULT;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.isPositionValid;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildIndelProbe;
import static com.hartwig.hmftools.panelbuilder.Utils.findDuplicates;

import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering a list of arbitrary small variants provided by the user.
public class CustomSmallVariants
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CUSTOM_SMALL_VARIANT;

    private static final Logger LOGGER = LogManager.getLogger(CustomSmallVariants.class);

    public static void generateProbes(final String variantFile, final RefGenomeInterface refGenome,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating custom small variant probes");

        List<CustomSmallVariant> customVariants = CustomSmallVariant.readFromFile(variantFile);

        checkVariantPositions(customVariants, refGenome.chromosomeLengths());
        checkRefPositionMatch(customVariants, refGenome);
        checkNoDuplicates(customVariants);

        generateProbes(customVariants, probeGenerator, panelData);

        LOGGER.info("Done generating custom small variant probes");
    }

    private static void checkVariantPositions(final List<CustomSmallVariant> customVariants, final Map<String, Integer> chromosomeLengths)
    {
        List<CustomSmallVariant> invalid = customVariants.stream()
                .filter(v -> !(isPositionValid(v.position(), chromosomeLengths)))
                .toList();
        if(!invalid.isEmpty())
        {
            invalid.forEach(v -> LOGGER.error("Invalid custom small variant positions: {}", v));
            throw new UserInputError("Invalid custom small variant positions");
        }
    }

    private static void checkRefPositionMatch(final List<CustomSmallVariant> customVariants, final RefGenomeInterface refGenome)
    {
        // Since the input format is redundant with the position and ref sequence, want to check that the ref sequence is actually correct.
        boolean valid = true;
        for(CustomSmallVariant v : customVariants)
        {
            String refGenomeSeq = refGenome.getBaseString(v.position().Chromosome, v.position().Position,
                    v.position().Position + v.refSequence().length() - 1);
            if(!v.refSequence().equals(refGenomeSeq))
            {
                LOGGER.error("Invalid custom small variant ref sequence. At {} expected {} but got {}", v.position(), refGenomeSeq, v.refSequence());
                valid = false;
            }
        }
        if(!valid)
        {
            throw new UserInputError("Invalid custom small variant ref sequence");
        }
    }

    private static void checkNoDuplicates(final List<CustomSmallVariant> customVariants)
    {
        LOGGER.debug("Checking custom small variants for duplicates");
        List<CustomSmallVariant> duplicated = findDuplicates(customVariants, (v1, v2) ->
                v1.position() == v2.position() && v1.refSequence().length() == v2.refSequence().length() &&
                        v1.altSequence().equals(v2.altSequence()));
        if(!duplicated.isEmpty())
        {
            duplicated.forEach(v -> LOGGER.error("Duplicate custom small variant: {}", v));
            throw new UserInputError("Duplicate custom small variants");
        }
    }

    private static void generateProbes(final List<CustomSmallVariant> customVariants, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        Stream<ProbeGenerationSpec> probeGenerationSpecs = customVariants.stream().map(CustomSmallVariants::createProbeGenerationSpec);
        probeGenerator.generateBatch(probeGenerationSpecs, panelData);
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final CustomSmallVariant customVariant)
    {
        LOGGER.debug("Generating probes for {}", customVariant);
        TargetMetadata metadata = new TargetMetadata(TARGET_TYPE, customVariant.extraInfo());
        SequenceDefinition definition =
                buildIndelProbe(customVariant.position().Chromosome, customVariant.position().Position, customVariant.refSequence(), customVariant.altSequence(), PROBE_LENGTH);
        // TODO: should target whole range or just altered bases?
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        ProbeEvaluator.Criteria evalCriteria = new ProbeEvaluator.Criteria(
                customVariant.qualityScoreMin() == null ? CUSTOM_SMALL_VARIANT_QUALITY_MIN_DEFAULT : customVariant.qualityScoreMin(),
                CUSTOM_SMALL_VARIANT_GC_TARGET, CUSTOM_SMALL_VARIANT_GC_TOLERANCE);
        return new ProbeGenerationSpec.SingleProbe(definition, targetedRange, metadata, evalCriteria);
    }
}
