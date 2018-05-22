package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ConsequenceDeterminer {

    private static final Logger LOGGER = LogManager.getLogger(ConsequenceDeterminer.class);

    @VisibleForTesting
    static final String FEATURE_TYPE_TRANSCRIPT = "transcript";

    @NotNull
    private final Set<String> transcripts;

    ConsequenceDeterminer(@NotNull final Set<String> transcripts) {
        this.transcripts = transcripts;
    }

    @NotNull
    List<VariantReport> run(@NotNull final List<SomaticVariant> variants) {
        Predicate<SomaticVariant> hasImpactingAnnotation = variant -> findImpactingAnnotation(variant, transcripts) != null;

        List<SomaticVariant> variantsToReport = variants.stream().filter(hasImpactingAnnotation).collect(Collectors.toList());

        return toVariantReport(variantsToReport);
    }

    @NotNull
    private List<VariantReport> toVariantReport(@NotNull final List<SomaticVariant> variantsToReport) {
        final List<VariantReport> reports = Lists.newArrayList();
        for (final SomaticVariant variant : variantsToReport) {
            final ImmutableVariantReport.Builder builder = ImmutableVariantReport.builder();
            final VariantAnnotation variantAnnotation = findImpactingAnnotation(variant, transcripts);
            // KODU: Variants with no impacting annotations should be filtered out by now.
            assert variantAnnotation != null;

            if (!variantAnnotation.allele().equals(variant.alt())) {
                LOGGER.warn("Annotated allele does not match alt from variant for " + variant);
            }

            builder.variant(variant);
            builder.gene(variantAnnotation.gene());
            builder.transcript(variantAnnotation.featureID());
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
    private static VariantAnnotation findImpactingAnnotation(@NotNull final SomaticVariant variant,
            @NotNull final Set<String> transcripts) {
        final List<VariantAnnotation> relevantAnnotations = findAllRelevantAnnotations(variant.annotations(), transcripts);
        for (final VariantAnnotation annotation : relevantAnnotations) {
            for (final VariantConsequence consequence : annotation.consequences()) {
                if (VariantConsequence.IMPACTING_CONSEQUENCES.contains(consequence)) {
                    return annotation;
                }
            }
        }
        return null;
    }

    @NotNull
    private static List<VariantAnnotation> findAllRelevantAnnotations(@NotNull final List<VariantAnnotation> annotations,
            @NotNull final Set<String> transcripts) {
        return annotations.stream()
                .filter(annotation -> annotation.featureType().equals(FEATURE_TYPE_TRANSCRIPT)
                        && transcripts.contains(annotation.featureID()))
                .collect(Collectors.toList());
    }
}
