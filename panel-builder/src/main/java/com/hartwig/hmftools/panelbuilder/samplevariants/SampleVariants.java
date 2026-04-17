package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantSelector.getDriverVariantCandidates;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantSelector.getNondriverVariantCandidates;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantSelector.markUnevaluatedCandidates;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantSelector.selectDriverVariants;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantSelector.selectNondriverVariants;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.panelbuilder.EvaluationResult;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.ProbeEvaluator;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationSpec;
import com.hartwig.hmftools.panelbuilder.ProbeGenerator;
import com.hartwig.hmftools.panelbuilder.SequenceDefinition;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.TargetedRange;
import com.hartwig.hmftools.panelbuilder.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Probes covering variants found in sample data.
// Inputs Linx and Purple data from a previous pipeline run.
// Probes are generated for driver SNV, INDEL, and SV, and then nondriver SNV/INDEL.
// Each variant gets 1 probe which consists of the alt sequence.
public class SampleVariants
{
    private static final ProbeEvaluator.Criteria NONDRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_NONDRIVER_QUALITY_MIN, SAMPLE_NONDRIVER_GC_TARGET, SAMPLE_NONDRIVER_GC_TOLERANCE);
    private static final ProbeEvaluator.Criteria DRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_DRIVER_QUALITY_MIN, SAMPLE_DRIVER_GC_TARGET, SAMPLE_DRIVER_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public record ExtraOutput(
            List<VariantInfo> variantInfos
    )
    {
    }

    public record VariantInfo(
            String variant,
            TargetMetadata.Type targetType,
            @Nullable VariantFilter filterReason
    )
    {
    }

    public static ExtraOutput generateProbes(final SampleVariantsConfig config, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating sample variant probes");

        List<Variant> variants = new ArrayList<>();

        checkSampleDirectories(config.purpleDir(), config.linxDir(), config.linxGermlineDir());
        variants.addAll(SomaticMutation.load(config.sampleId(), config.purpleDir()));
        variants.addAll(GermlineMutation.load(config.sampleId(), config.purpleDir()));
        variants.addAll(SomaticSv.load(config.sampleId(), config.purpleDir(), config.linxDir()));
        variants.addAll(GermlineSv.load(config.sampleId(), config.linxGermlineDir()));

        Map<Variant, EvaluationResult<VariantFilter>> variantFilters = new HashMap<>();
        generateProbes(variants, config.maxProbes(), config.prioritiseSmallIndels(), probeGenerator, panelData, variantFilters);

        List<VariantInfo> variantInfos = createVariantInfos(variants, variantFilters);
        ExtraOutput extraOutput = new ExtraOutput(variantInfos);

        LOGGER.info("Done generating sample variant probes");

        return extraOutput;
    }

    private static void checkSampleDirectories(final String purpleDir, final String linxDir, final String linxGermlineDir)
    {
        // Allow Linx inputs to be optional.
        if(!Files.exists(Paths.get(purpleDir))
                || (linxDir != null && !Files.exists(Paths.get(linxDir)))
                || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            throw new UserInputError("Missing Purple or Linx directories");
        }
    }

    public static void generateProbes(final List<Variant> variants, int maxProbes, boolean prioritiseSmallIndels,
            final ProbeGenerator probeGenerator, PanelData panelData, Map<Variant, EvaluationResult<VariantFilter>> variantFilters)
    {
        // Variants which can be considered for probe generation.
        List<Variant> driverVariantCandidates = getDriverVariantCandidates(variants);
        List<SomaticMutation> nonDriverVariantCandidates = getNondriverVariantCandidates(variants, prioritiseSmallIndels);

        int generatedProbes = 0;
        Map<String, Integer> geneDisruptions = new HashMap<>();
        generatedProbes += generateProbes(
                remainingProbes -> selectDriverVariants(driverVariantCandidates, geneDisruptions, variantFilters, remainingProbes),
                maxProbes - generatedProbes, probeGenerator, panelData);
        generatedProbes += generateProbes(
                remainingProbes -> selectNondriverVariants(nonDriverVariantCandidates, variantFilters, remainingProbes),
                maxProbes - generatedProbes, probeGenerator, panelData);

        // At this point, candidates are either evaluated (accepted or rejected) or couldn't be evaluated due to the max probe count.
        markUnevaluatedCandidates(driverVariantCandidates, variantFilters);
        markUnevaluatedCandidates(nonDriverVariantCandidates, variantFilters);
    }

    private static int generateProbes(final Function<Integer, List<Variant>> variantSelector, int maxProbes,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        // Want to create exactly maxProbes probes if there are enough variants, but since probe evaluation is done in batches and can
        // reject probes, we need to iteratively generate until the quota is filled.

        int generatedProbes = 0;
        while(true)
        {
            int remainingProbes = maxProbes - generatedProbes;
            if(remainingProbes <= 0)
            {
                // Filled the probe quota.
                break;
            }

            List<Variant> variants = variantSelector.apply(remainingProbes);
            if(variants.isEmpty())
            {
                // Can generate more probes, but there are no more candidate variants available to select.
                break;
            }

            Stream<ProbeGenerationSpec> probeGenerationSpecs = variants.stream().map(SampleVariants::createProbeGenerationSpec);
            generatedProbes += probeGenerator.generateBatch(probeGenerationSpecs, panelData).probes().size();
        }
        return generatedProbes;
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final Variant variant)
    {
        LOGGER.trace("Generating probe for variant: {}", variant);

        SequenceDefinition definition = variant.generateProbe();
        // TODO: should target whole range or just altered bases?
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        TargetMetadata metadata = createTargetMetadata(variant);
        ProbeEvaluator.Criteria evalCriteria = variant.isDriver() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
        return new ProbeGenerationSpec.SingleProbe(definition, targetedRange, metadata, evalCriteria);
    }

    private static TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(variant.targetType(), variant.toString());
    }

    private static List<VariantInfo> createVariantInfos(final List<Variant> variants,
            final Map<Variant, EvaluationResult<VariantFilter>> variantFilters)
    {
        return variants.stream()
                .map(variant ->
                {
                    EvaluationResult<VariantFilter> filterResult = variantFilters.get(variant);
                    // If we didn't even attempt to filter the variant, it must have been because it was never a candidate.
                    VariantFilter filterReason = filterResult == null ? VariantFilter.NotCandidate : filterResult.rejectionInfo();
                    return new VariantInfo(variant.toString(), variant.targetType(), filterReason);
                })
                .toList();
    }
}
