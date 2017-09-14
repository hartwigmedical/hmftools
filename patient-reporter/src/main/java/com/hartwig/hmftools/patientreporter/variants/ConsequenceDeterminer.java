package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.and;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ConsequenceDeterminer {

    private static final Logger LOGGER = LogManager.getLogger(ConsequenceDeterminer.class);

    @VisibleForTesting
    static final String FEATURE_TYPE_TRANSCRIPT = "transcript";

    @NotNull
    private final Slicer hmfSlicingRegion;
    @NotNull
    private final Map<String, HmfGenomeRegion> relevantTranscriptMap;

    ConsequenceDeterminer(@NotNull final GeneModel geneModel) {
        this(geneModel.slicer(), geneModel.transcriptMap());
    }

    ConsequenceDeterminer(@NotNull final Slicer hmfSlicingRegion, @NotNull final Map<String, HmfGenomeRegion> relevantTranscriptMap) {
        this.hmfSlicingRegion = hmfSlicingRegion;
        this.relevantTranscriptMap = relevantTranscriptMap;
    }

    @NotNull
    ConsequenceOutput run(@NotNull final List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> consequenceRule =
                and(isIncludedIn(hmfSlicingRegion), hasActionableConsequence(relevantTranscriptMap.keySet()));

        final List<SomaticVariant> consequentialVariants = filter(variants, consequenceRule);

        return new ConsequenceOutput(consequentialVariants, toVariantReport(consequentialVariants));
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull final Slicer slicer) {
        return slicer::includes;
    }

    @NotNull
    private static Predicate<SomaticVariant> hasActionableConsequence(@NotNull final Set<String> relevantTranscripts) {
        return variant -> findPrimaryRelevantAnnotation(variant, relevantTranscripts) != null;
    }

    @NotNull
    private List<VariantReport> toVariantReport(@NotNull final List<SomaticVariant> variants) {
        final List<VariantReport> reports = Lists.newArrayList();
        for (final SomaticVariant variant : variants) {
            final ImmutableVariantReport.Builder builder = ImmutableVariantReport.builder();
            final VariantAnnotation variantAnnotation = findPrimaryRelevantAnnotation(variant, relevantTranscriptMap.keySet());
            // KODU: Variants with no relevant annotations should be filtered out by now.
            assert variantAnnotation != null;

            final HmfGenomeRegion hmfGenomeRegion = relevantTranscriptMap.get(variantAnnotation.featureID());
            assert hmfGenomeRegion != null;

            if (!variantAnnotation.gene().equals(hmfGenomeRegion.gene())) {
                LOGGER.warn("Annotated gene does not match gene expected from slicing annotation for " + variant);
            }

            if (!variantAnnotation.allele().equals(variant.alt())) {
                LOGGER.warn("Annotated allele does not match alt from variant for " + variant);
            }

            builder.gene(variantAnnotation.gene());
            builder.chromosome(variant.chromosome());
            builder.position(variant.position());
            builder.ref(variant.ref());
            builder.alt(variant.alt());
            builder.transcript(hmfGenomeRegion.transcript());
            builder.hgvsCoding(variantAnnotation.hgvsCoding());
            builder.hgvsProtein(variantAnnotation.hgvsProtein());
            builder.consequence(variantAnnotation.consequenceString());
            final String cosmicID = variant.cosmicID();
            if (cosmicID != null) {
                builder.cosmicID(cosmicID);
            }
            builder.totalReadCount(variant.totalReadCount());
            builder.alleleReadCount(variant.alleleReadCount());
            reports.add(builder.build());
        }

        return reports;
    }

    @Nullable
    private static VariantAnnotation findPrimaryRelevantAnnotation(@NotNull final SomaticVariant variant,
            @NotNull final Set<String> relevantTranscripts) {
        final List<VariantAnnotation> relevantAnnotations = findAllRelevantAnnotations(variant.annotations(), relevantTranscripts);
        for (final VariantAnnotation annotation : relevantAnnotations) {
            for (final VariantConsequence consequence : annotation.consequences())
                if (VariantConsequence.ACTIONABLE_CONSEQUENCES.contains(consequence)) {
                    return annotation;
                }
        }
        return null;
    }

    @NotNull
    private static List<VariantAnnotation> findAllRelevantAnnotations(@NotNull final List<VariantAnnotation> annotations,
            @NotNull final Set<String> relevantTranscripts) {
        return annotations.stream()
                .filter(annotation -> annotation.featureType().equals(FEATURE_TYPE_TRANSCRIPT) && relevantTranscripts.contains(
                        annotation.featureID()))
                .collect(Collectors.toList());
    }
}
