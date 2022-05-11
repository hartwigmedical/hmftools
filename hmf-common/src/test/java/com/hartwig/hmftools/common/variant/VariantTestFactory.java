package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VariantTestFactory {

    private VariantTestFactory() {
    }

    @NotNull
    public static ReportableVariant create() {
        return builder().build();
    }

    @NotNull
    public static ImmutableReportableVariant.Builder builder() {
        return ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene(Strings.EMPTY)
                .transcript("transcript")
                .isCanonical(true)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.SNP)
                .otherReportedEffects(Strings.EMPTY)
                .canonicalTranscript("123")
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .minorAlleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false);
    }
}
