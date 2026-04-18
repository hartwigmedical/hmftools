package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INDEL_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INSERT_SEQUENCE_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SV_BREAKENDS_PER_GENE_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.panelbuilder.EvaluationResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VariantSelector
{
    private static final Logger LOGGER = LogManager.getLogger(VariantSelector.class);

    public static List<Variant> selectDriverVariants(final List<Variant> candidateVariants, Map<String, Integer> geneDisruptions,
            Map<Variant, EvaluationResult<VariantFilter>> variantFilters, int maxCount)
    {
        LOGGER.debug("Selecting up to {} driver variants", maxCount);

        List<Variant> selectedVariants = new ArrayList<>();
        for(Variant variant : candidateVariants)
        {
            if(hasBeenEvaluated(variant, variantFilters))
            {
                continue;
            }

            if(selectedVariants.size() >= maxCount)
            {
                // We assume that there will be enough probes to cover all the driver variants.
                // If not, then we randomly(?) discarded some drivers, which may be important to handle.
                // Potentially, there should be a variant prioritisation scheme to avoid this scenario.
                LOGGER.warn("Filled sample variant probe quota without including all driver variants!");
                break;
            }

            EvaluationResult<VariantFilter> filterResult = driverFilters(variant, geneDisruptions);
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

    public static List<Variant> getDriverVariantCandidates(final List<Variant> variants)
    {
        List<Variant> candidates = variants.stream()
                .filter(Variant::isDriver)
                .toList();
        LOGGER.debug("{} driver variant candidates out of {} variants", candidates.size(), variants.size());
        return candidates;
    }

    public static List<Variant> selectNondriverVariants(final List<SomaticMutation> candidateVariants,
            Map<Variant, EvaluationResult<VariantFilter>> variantFilters, int maxCount)
    {
        LOGGER.debug("Selecting up to {} nondriver variants", maxCount);

        List<Variant> selectedVariants = new ArrayList<>();
        for(SomaticMutation variant : candidateVariants)
        {
            if(selectedVariants.size() >= maxCount)
            {
                break;
            }

            if(hasBeenEvaluated(variant, variantFilters))
            {
                continue;
            }

            EvaluationResult<VariantFilter> filterResult = nondriverFilters(variant);
            variantFilters.put(variant, filterResult);
            filterResult.unwrap(
                    () -> selectedVariants.add(variant),
                    failReason -> LOGGER.trace("Variant failed filter: {{}} reason=\"{}\"", variant, failReason)
            );
        }
        return selectedVariants;
    }

    public static List<SomaticMutation> getNondriverVariantCandidates(final List<Variant> variants, boolean prioritiseSmallIndels)
    {
        // For nondrivers, we are only interested in somatic SNV/INDEL.
        // Also, some variants are prioritised over others.
        List<SomaticMutation> candidates = variants.stream()
                .filter(variant -> variant instanceof SomaticMutation)
                .filter(variant -> !variant.isDriver())
                .map(variant -> (SomaticMutation) variant)
                .sorted(new NondriverVariantComparator(prioritiseSmallIndels))
                .toList();
        LOGGER.debug("{} nondriver variant candidates out of {} variants", candidates.size(), variants.size());
        return candidates;
    }

    private record NondriverVariantComparator(
            boolean prioritiseSmallerVariants
    )
            implements Comparator<SomaticMutation>
    {
        @Override
        public int compare(final SomaticMutation v1, final SomaticMutation v2)
        {
            if(prioritiseSmallerVariants && v1.indelLength() != v2.indelLength())
            {
                return v1.indelLength() < v2.indelLength() ? -1 : 1;
            }
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

    private static boolean hasBeenEvaluated(final Variant variant, final Map<Variant, EvaluationResult<VariantFilter>> variantFilters)
    {
        return variantFilters.containsKey(variant);
    }

    private static EvaluationResult<VariantFilter> commonFilters(final Variant variant)
    {
        List<Supplier<EvaluationResult<VariantFilter>>> filters = new ArrayList<>();

        if(variant instanceof StructuralVariant sv)
        {
            filters.add(() -> insertSequenceLengthFilter(sv.insertSequenceLength()));
        }

        return EvaluationResult.applyEvaluations(filters);
    }

    private static EvaluationResult<VariantFilter> driverFilters(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        return EvaluationResult.applyEvaluations(List.of(
                () -> commonFilters(variant),
                () -> driverOnlyFilters(variant, geneDisruptions)));
    }

    private static EvaluationResult<VariantFilter> driverOnlyFilters(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        List<Supplier<EvaluationResult<VariantFilter>>> filters = new ArrayList<>();

        if(variant instanceof SomaticSv sv)
        {
            if(sv.isReportedDisruption())
            {
                filters.add(() -> vafFilter(sv.vaf()));
                filters.add(() -> tumorFragmentsFilter(sv.tumorFragments()));
            }
        }

        filters.add(() -> geneDisruptionFilter(variant, geneDisruptions));

        return EvaluationResult.applyEvaluations(filters);
    }

    private static EvaluationResult<VariantFilter> nondriverFilters(final SomaticMutation variant)
    {
        return EvaluationResult.applyEvaluations(List.of(
                () -> commonFilters(variant),
                () -> nondriverOnlyFilters(variant)));
    }

    private static EvaluationResult<VariantFilter> nondriverOnlyFilters(final SomaticMutation variant)
    {
        return EvaluationResult.applyEvaluations(List.of(
                () -> vafFilter(variant.vaf()),
                () -> tumorFragmentsFilter(variant.tumorFragments()),
                () -> indelLengthFilter(variant.indelLength()),
                () -> repeatCountFilter(variant.repeatCount()),
                () -> germlineStatusFilter(variant.germlineStatus())));
    }

    private static EvaluationResult<VariantFilter> insertSequenceLengthFilter(int insertSequenceLength)
    {
        return EvaluationResult.condition(insertSequenceLength <= SAMPLE_INSERT_SEQUENCE_LENGTH_MAX, VariantFilter.InsertSeqLength);
    }

    private static EvaluationResult<VariantFilter> vafFilter(double vaf)
    {
        return EvaluationResult.condition(vaf >= SAMPLE_VAF_MIN, VariantFilter.VAF);
    }

    private static EvaluationResult<VariantFilter> tumorFragmentsFilter(int tumorFragments)
    {
        return EvaluationResult.condition(tumorFragments >= SAMPLE_FRAGMENT_COUNT_MIN, VariantFilter.TumorFragments);
    }

    private static EvaluationResult<VariantFilter> indelLengthFilter(int indelLength)
    {
        return EvaluationResult.condition(indelLength <= SAMPLE_INDEL_LENGTH_MAX, VariantFilter.IndelLength);
    }

    private static EvaluationResult<VariantFilter> repeatCountFilter(int repeatCount)
    {
        return EvaluationResult.condition(repeatCount <= SAMPLE_REPEAT_COUNT_MAX, VariantFilter.RepeatCount);
    }

    private static EvaluationResult<VariantFilter> germlineStatusFilter(final GermlineStatus germlineStatus)
    {
        return EvaluationResult.condition(germlineStatus == GermlineStatus.DIPLOID, VariantFilter.GermlineStatus);
    }

    private static EvaluationResult<VariantFilter> geneDisruptionFilter(final Variant variant,
            final Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            boolean pass = sv.disruptedGenes().stream()
                    .allMatch(gene -> geneDisruptions.getOrDefault(gene, 0) + 1 <= SAMPLE_SV_BREAKENDS_PER_GENE_MAX);
            return EvaluationResult.condition(pass, VariantFilter.GeneDisruption);
        }
        return EvaluationResult.accept();
    }

    private static void registerDisruptedGenes(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            sv.disruptedGenes().forEach(gene -> geneDisruptions.put(gene, geneDisruptions.getOrDefault(gene, 0) + 1));
        }
    }

    public static <T extends Variant> void markUnevaluatedCandidates(final List<T> candidates,
            Map<Variant, EvaluationResult<VariantFilter>> filters)
    {
        List<T> unevaluated = candidates.stream()
                .filter(candidate -> !hasBeenEvaluated(candidate, filters))
                .toList();
        for(Variant variant : unevaluated)
        {
            filters.put(variant, EvaluationResult.reject(VariantFilter.ProbeLimit));
        }
    }
}
