package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
    List<VariantReport> run(@NotNull final List<EnrichedSomaticVariant> variants) {
        Predicate<EnrichedSomaticVariant> hasImpactingAnnotation = variant -> findImpactingAnnotation(variant, transcripts) != null;

        List<EnrichedSomaticVariant> variantsToReport = variants.stream().filter(hasImpactingAnnotation).collect(Collectors.toList());

        return toVariantReport(variantsToReport);
    }

    @NotNull
    private List<VariantReport> toVariantReport(@NotNull final List<EnrichedSomaticVariant> variantsToReport) {
        final List<VariantReport> reports = Lists.newArrayList();
        for (final EnrichedSomaticVariant variant : variantsToReport) {
            final ImmutableVariantReport.Builder builder = ImmutableVariantReport.builder();
            final SnpEffAnnotation snpEffAnnotation = findImpactingAnnotation(variant, transcripts);
            // KODU: Variants with no impacting annotations should be filtered out by now.
            assert snpEffAnnotation != null;

            if (!snpEffAnnotation.allele().equals(variant.alt())) {
                LOGGER.warn("Annotated allele does not match alt from variant for " + variant);
            }

            builder.variant(variant);
            builder.gene(snpEffAnnotation.gene());
            builder.proteinImpact(snpEffAnnotation.hgvsProtein());
            builder.proteinImpactType(snpEffAnnotation.consequenceString());
            final String cosmicID = variant.canonicalCosmicID();
            if (cosmicID != null) {
                builder.knowledgebaseKey(cosmicID);
            }
            builder.totalReadCount(variant.totalReadCount());
            builder.alleleReadCount(variant.alleleReadCount());
            builder.ploidy(Strings.EMPTY);
            builder.purityAdjustedVAF(0D);
            builder.clonalityStatus(Strings.EMPTY);
            builder.wildTypeStatus(Strings.EMPTY);
            builder.driverStatus(Strings.EMPTY);
            builder.actionabilityStatus(Strings.EMPTY);
            reports.add(builder.build());
        }

        return reports;
    }

    @Nullable
    private static SnpEffAnnotation findImpactingAnnotation(@NotNull final SomaticVariant variant, @NotNull final Set<String> transcripts) {
        final List<SnpEffAnnotation> relevantAnnotations = findAllRelevantAnnotations(variant.snpEffAnnotations(), transcripts);
        for (final SnpEffAnnotation annotation : relevantAnnotations) {
            for (final VariantConsequence consequence : annotation.consequences()) {
                if (VariantConsequence.IMPACTING_CONSEQUENCES.contains(consequence)) {
                    return annotation;
                }
            }
        }
        return null;
    }

    @NotNull
    private static List<SnpEffAnnotation> findAllRelevantAnnotations(@NotNull final List<SnpEffAnnotation> annotations,
            @NotNull final Set<String> transcripts) {
        return annotations.stream()
                .filter(annotation -> annotation.featureType().equals(FEATURE_TYPE_TRANSCRIPT)
                        && transcripts.contains(annotation.transcript()))
                .collect(Collectors.toList());
    }
}
