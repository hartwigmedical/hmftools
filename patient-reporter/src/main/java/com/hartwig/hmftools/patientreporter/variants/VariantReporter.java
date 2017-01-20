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
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class VariantReporter {

    private static final Logger LOGGER = LogManager.getLogger(VariantReporter.class);

    @VisibleForTesting
    static final String FEATURE_TYPE_TRANSCRIPT = "transcript";
    // KODU: This boolean exists to evaluate the impact of annotation-filtering on actual patients.
    private static final boolean INCLUDE_ALL_ANNOTATIONS_FOR_IMPACT = false;

    private static final List<VariantConsequence> ACTIONABLE_CONSEQUENCES = Lists.newArrayList(
            VariantConsequence.TRANSCRIPT_ABLATION, VariantConsequence.TRANSCRIPT_AMPLIFICATION,
            VariantConsequence.SPLICE_ACCEPTOR_VARIANT, VariantConsequence.SPLICE_DONOR_VARIANT,
            VariantConsequence.SPLICE_REGION_VARIANT, VariantConsequence.STOP_GAINED, VariantConsequence.STOP_LOST,
            VariantConsequence.INCOMPLETE_TERMINAL_CODING_VARIANT, VariantConsequence.INITIATOR_CODON_VARIANT,
            VariantConsequence.START_LOST, VariantConsequence.FRAMESHIFT_VARIANT, VariantConsequence.INFRAME_INSERTION,
            VariantConsequence.INFRAME_DELETION, VariantConsequence.MISSENSE_VARIANT);

    @NotNull
    private final Slicer hmfSlicingRegion;
    @NotNull
    private final Map<String, HMFSlicingAnnotation> relevantTranscriptMap;

    VariantReporter(@NotNull final Slicer hmfSlicingRegion,
            @NotNull final Map<String, HMFSlicingAnnotation> relevantTranscriptMap) {
        this.hmfSlicingRegion = hmfSlicingRegion;
        this.relevantTranscriptMap = relevantTranscriptMap;
    }

    @NotNull
    List<VariantReport> run(@NotNull final List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> consequenceRule = and(isIncludedIn(hmfSlicingRegion),
                hasActionableConsequence(relevantTranscriptMap.keySet()));

        final List<SomaticVariant> variantsWithConsequence = filter(variants, consequenceRule);

        return toVariantReport(variantsWithConsequence);
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
            final VariantReport.Builder builder = new VariantReport.Builder();
            final VariantAnnotation variantAnnotation = findPrimaryRelevantAnnotation(variant,
                    relevantTranscriptMap.keySet());
            // KODU: Variants with no relevant annotations should be filtered out by now.
            assert variantAnnotation != null;

            final HMFSlicingAnnotation slicingAnnotation = relevantTranscriptMap.get(variantAnnotation.featureID());
            assert slicingAnnotation != null;

            if (!variantAnnotation.gene().equals(slicingAnnotation.gene())) {
                LOGGER.warn("Annotated gene does not match gene expected from slicing annotation for " + variant);
            }

            if (!variantAnnotation.allele().equals(variant.alt())) {
                LOGGER.warn("Annotated allele does not match alt from variant for " + variant);
            }

            builder.gene(variantAnnotation.gene());
            builder.position(variant.chromosome() + ":" + variant.position());
            builder.ref(variant.ref());
            builder.alt(variant.alt());
            builder.transcript(slicingAnnotation.transcriptID() + "." + slicingAnnotation.transcriptVersion());
            builder.hgvsCoding(variantAnnotation.hgvsCoding());
            builder.hgvsProtein(variantAnnotation.hgvsProtein());
            builder.consequence(variantAnnotation.consequenceString());
            final String cosmicID = variant.cosmicID();
            if (cosmicID != null) {
                builder.cosmicID(cosmicID);
            }
            builder.alleleFrequency(Double.toString(variant.alleleFrequency()));
            builder.readDepth(Integer.toString(variant.readDepth()));
            reports.add(builder.build());
        }

        return reports;
    }

    @Nullable
    private static VariantAnnotation findPrimaryRelevantAnnotation(@NotNull final SomaticVariant variant,
            @NotNull final Set<String> relevantTranscripts) {
        final List<VariantAnnotation> relevantAnnotations = findAllRelevantAnnotations(variant.annotations(),
                relevantTranscripts);
        for (final VariantAnnotation annotation : relevantAnnotations) {
            for (final VariantConsequence consequence : annotation.consequences())
                if (ACTIONABLE_CONSEQUENCES.contains(consequence)) {
                    return annotation;
                }
        }
        return null;
    }

    @NotNull
    private static List<VariantAnnotation> findAllRelevantAnnotations(
            @NotNull final List<VariantAnnotation> annotations, @NotNull final Set<String> relevantTranscripts) {
        if (INCLUDE_ALL_ANNOTATIONS_FOR_IMPACT) {
            return annotations;
        }

        return annotations.stream().filter(
                annotation -> annotation.featureType().equals(FEATURE_TYPE_TRANSCRIPT) && relevantTranscripts.contains(
                        annotation.featureID())).collect(Collectors.toList());
    }
}
