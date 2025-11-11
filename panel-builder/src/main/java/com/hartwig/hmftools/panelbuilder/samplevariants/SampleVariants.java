package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INDEL_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SV_BREAKENDS_PER_GENE_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.ProbeEvaluator;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
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
            @Nullable String filterReason
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

        Map<Variant, FilterResult> variantFilters = new HashMap<>();
        ProbeGenerationResult result = generateProbes(variants, config.maxProbes(), probeGenerator, panelData, variantFilters);

        List<VariantInfo> variantInfos = createVariantInfos(variants, variantFilters);
        ExtraOutput extraOutput = new ExtraOutput(variantInfos);

        panelData.addResult(result);

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

    public static ProbeGenerationResult generateProbes(final List<Variant> variants, int maxProbes, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage, Map<Variant, FilterResult> variantFilters)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();
        Map<String, Integer> geneDisruptions = new HashMap<>();
        result = result.add(
                generateProbes(
                        remainingProbes -> selectDriverVariants(variants, geneDisruptions, variantFilters, remainingProbes),
                        maxProbes - result.probes().size(), probeGenerator, coverage));
        result = result.add(
                generateProbes(
                        remainingProbes -> selectNondriverVariants(variants, variantFilters, remainingProbes),
                        maxProbes - result.probes().size(), probeGenerator, coverage));
        return result;
    }

    private static ProbeGenerationResult generateProbes(final Function<Integer, List<Variant>> variantSelector, int maxProbes,
            final ProbeGenerator probeGenerator, final PanelCoverage coverage)
    {
        // Want to create exactly maxProbes probes if there are enough variants, but since probe evaluation is done in batches and can
        // reject probes, we need to iteratively generate until the quota is filled.

        ProbeGenerationResult result = new ProbeGenerationResult();
        while(true)
        {
            int remainingProbes = maxProbes - result.probes().size();
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
            result = result.add(probeGenerator.generateBatch(probeGenerationSpecs, coverage));
        }
        return result;
    }

    private static List<Variant> selectDriverVariants(final List<Variant> variants, Map<String, Integer> geneDisruptions,
            Map<Variant, FilterResult> variantFilters, int maxCount)
    {
        LOGGER.debug("Selecting up to {} driver variants", maxCount);

        List<Variant> candidateVariants = variants.stream()
                .filter(variant -> !variantFilters.containsKey(variant))
                .filter(Variant::isDriver)
                .toList();

        List<Variant> selectedVariants = new ArrayList<>();
        for(Variant variant : candidateVariants)
        {
            if(selectedVariants.size() >= maxCount)
            {
                // We assume that there will be enough probes to cover all the driver variants.
                // If not, then we randomly(?) discarded some drivers, which may be important to handle.
                // Potentially, there should be a variant prioritisation scheme to avoid this scenario.
                LOGGER.warn("Filled sample variant probe quota without including all driver variants!");
                break;
            }

            FilterResult filterResult = driverFilters(variant, geneDisruptions);
            variantFilters.put(variant, filterResult);
            filterResult.unwrap(
                    () ->
                    {
                        registerDisruptedGenes(variant, geneDisruptions);
                        selectedVariants.add(variant);
                    },
                    failReason -> LOGGER.trace("Variant failed filter: {{}} reason=\"{}\"", variant, failReason)
            );
        }
        return selectedVariants;
    }

    private static List<Variant> selectNondriverVariants(final List<Variant> variants, Map<Variant, FilterResult> variantFilters,
            int maxCount)
    {
        LOGGER.debug("Selecting up to {} nondriver variants", maxCount);

        // For nondrivers, we are only interested in somatic SNV/INDEL.
        // Also, some variants are prioritised over others.
        List<SomaticMutation> candidateVariants = variants.stream()
                .filter(variant -> !variantFilters.containsKey(variant))
                .filter(variant -> variant instanceof SomaticMutation)
                .filter(variant -> !variant.isDriver())
                .map(variant -> (SomaticMutation) variant)
                .sorted(new NondriverVariantComparator())
                .toList();

        List<Variant> selectedVariants = new ArrayList<>();
        for(SomaticMutation variant : candidateVariants)
        {
            if(selectedVariants.size() >= maxCount)
            {
                break;
            }

            FilterResult filterResult = nondriverFilters(variant);
            variantFilters.put(variant, filterResult);
            filterResult.unwrap(
                    () -> selectedVariants.add(variant),
                    failReason -> LOGGER.trace("Variant failed filter: {{}} reason=\"{}\"", variant, failReason)
            );
        }
        return selectedVariants;
    }

    private static FilterResult driverFilters(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        List<Supplier<FilterResult>> filters = new ArrayList<>();

        if(variant instanceof SomaticSv sv)
        {
            if(sv.isReportedDisruption())
            {
                filters.add(() -> vafFilter(sv.vaf()));
                filters.add(() -> tumorFragmentsFilter(sv.tumorFragments()));
            }
        }

        filters.add(() -> geneDisruptionFilter(variant, geneDisruptions));

        return FilterResult.applyFilters(filters);
    }

    private static class NondriverVariantComparator implements Comparator<SomaticMutation>
    {
        @Override
        public int compare(final SomaticMutation v1, final SomaticMutation v2)
        {
            if(v1.isCoding() != v2.isCoding())
            {
                return v1.isCoding() ? -1 : 1;
            }
            if(v1.isClonal() != v2.isClonal())
            {
                return v1.isClonal() ? -1 : 1;
            }
            // Otherwise random but deterministic ordering.
            return Integer.compare(v1.deterministicHash(), v2.deterministicHash());
        }
    }

    private static FilterResult nondriverFilters(final SomaticMutation variant)
    {
        return FilterResult.applyFilters(List.of(
                () -> vafFilter(variant.vaf()),
                () -> tumorFragmentsFilter(variant.tumorFragments()),
                () -> indelLengthFilter(variant.indelLength()),
                () -> repeatCountFilter(variant.repeatCount()),
                () -> germlineStatusFilter(variant.germlineStatus())));
    }

    private static FilterResult vafFilter(double vaf)
    {
        return FilterResult.condition(vaf >= SAMPLE_VAF_MIN, "VAF");
    }

    private static FilterResult tumorFragmentsFilter(int tumorFragments)
    {
        return FilterResult.condition(tumorFragments >= SAMPLE_FRAGMENT_COUNT_MIN, "tumor fragments");
    }

    private static FilterResult indelLengthFilter(int indelLength)
    {
        return FilterResult.condition(indelLength <= SAMPLE_INDEL_LENGTH_MAX, "indel length");
    }

    private static FilterResult repeatCountFilter(int repeatCount)
    {
        return FilterResult.condition(repeatCount <= SAMPLE_REPEAT_COUNT_MAX, "repeat count");
    }

    private static FilterResult germlineStatusFilter(final GermlineStatus germlineStatus)
    {
        return FilterResult.condition(germlineStatus == GermlineStatus.DIPLOID, "germline status");
    }

    private static FilterResult geneDisruptionFilter(final Variant variant, final Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            boolean pass = sv.disruptedGenes().stream()
                    .allMatch(gene -> geneDisruptions.getOrDefault(gene, 0) + 1 <= SAMPLE_SV_BREAKENDS_PER_GENE_MAX);
            return FilterResult.condition(pass, "gene disruptions");
        }
        return FilterResult.pass();
    }

    private static void registerDisruptedGenes(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            sv.disruptedGenes().forEach(gene -> geneDisruptions.put(gene, geneDisruptions.getOrDefault(gene, 0) + 1));
        }
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final Variant variant)
    {
        LOGGER.trace("Generating probe for variant: {}", variant);

        SequenceDefinition definition = variant.generateProbe();
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        TargetMetadata metadata = createTargetMetadata(variant);
        ProbeEvaluator.Criteria evalCriteria = variant.isDriver() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
        return new ProbeGenerationSpec.SingleProbe(definition, targetedRange, metadata, evalCriteria);
    }

    private static TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(variant.targetType(), variant.toString());
    }

    private static List<VariantInfo> createVariantInfos(final List<Variant> variants, final Map<Variant, FilterResult> variantFilters)
    {
        return variants.stream()
                .map(variant ->
                {
                    FilterResult filterResult = variantFilters.get(variant);
                    // If we didn't even attempt to filter the variant, it must have been because we encountered the probe limit first.
                    String filterReason = filterResult == null ? "max probe count" : filterResult.failReason();
                    return new VariantInfo(variant.toString(), filterReason);
                })
                .toList();
    }
}
