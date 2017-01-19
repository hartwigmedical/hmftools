package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.and;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class ConsequenceRule {

    private static final Logger LOGGER = LogManager.getLogger(ConsequenceRule.class);

    private static final String FEATURE_TYPE_TRANSCRIPT = "transcript";
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
    private final List<String> relevantTranscripts;

    @NotNull
    static ConsequenceRule fromHmfSlicingRegion(@NotNull final Slicer hmfSlicingRegion) {
        return new ConsequenceRule(hmfSlicingRegion, extractTranscripts(hmfSlicingRegion));
    }

    @NotNull
    private static List<String> extractTranscripts(final @NotNull Slicer hmfSlicingRegion) {
        final List<String> transcripts = Lists.newArrayList();
        for (final GenomeRegion region : hmfSlicingRegion.regions()) {
            HMFSlicingAnnotation annotation = HMFSlicingAnnotation.fromGenomeRegion(region);
            if (annotation == null) {
                LOGGER.warn("Could not extract annotation from hmf slicing region: " + region);
            } else {
                transcripts.add(annotation.transcriptID());
            }
        }
        return transcripts;
    }

    private ConsequenceRule(@NotNull final Slicer hmfSlicingRegion, @NotNull final List<String> relevantTranscripts) {
        this.hmfSlicingRegion = hmfSlicingRegion;
        this.relevantTranscripts = relevantTranscripts;
    }

    @NotNull
    List<SomaticVariant> apply(@NotNull List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> consequenceRule = and(isIncludedIn(hmfSlicingRegion),
                hasActionableConsequence(relevantTranscripts));

        return filter(variants, consequenceRule);
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }

    @NotNull
    private static Predicate<SomaticVariant> hasActionableConsequence(
            @NotNull final List<String> relevantTranscripts) {
        return variant -> {
            final List<VariantAnnotation> annotations = findRelevantAnnotations(variant.annotations(),
                    relevantTranscripts);
            if (annotations.isEmpty()) {
                LOGGER.warn("Variants should be filtered on transcripts already by this stage!");
                return false;
            }
            for (final VariantAnnotation annotation : annotations) {
                for (VariantConsequence consequence : annotation.consequences()) {
                    if (ACTIONABLE_CONSEQUENCES.contains(consequence)) {
                        return true;
                    }
                }
            }
            return false;
        };
    }

    @NotNull
    private static List<VariantAnnotation> findRelevantAnnotations(final List<VariantAnnotation> annotations,
            final List<String> relevantTranscripts) {
        if (INCLUDE_ALL_ANNOTATIONS_FOR_IMPACT) {
            return annotations;
        }

        return annotations.stream().filter(
                annotation -> annotation.featureType().equals(FEATURE_TYPE_TRANSCRIPT) && relevantTranscripts.contains(
                        annotation.featureID())).collect(Collectors.toList());
    }
}
