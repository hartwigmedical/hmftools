package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;

import org.jetbrains.annotations.NotNull;

public class PurpleConversion {

    private PurpleConversion() {
    }

    @NotNull
    public static PurpleCopyNumber convert(com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .averageTumorCopyNumber(copyNumber.averageTumorCopyNumber())
                .build();
    }

    @NotNull
    public static PurpleGeneCopyNumber convert(GeneCopyNumber geneCopyNumber) {
        return ImmutablePurpleGeneCopyNumber.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .geneName(geneCopyNumber.geneName())
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .minMinorAlleleCopyNumber(geneCopyNumber.minMinorAlleleCopyNumber())
                .build();
    }

    @NotNull
    public static PurpleDriver convert(DriverCatalog catalog) {
        return ImmutablePurpleDriver.builder()
                .gene(catalog.gene())
                .transcript(catalog.transcript())
                .driver(PurpleDriverType.valueOf(catalog.driver().name()))
                .driverLikelihood(catalog.driverLikelihood())
                .likelihoodMethod(PurpleLikelihoodMethod.valueOf(catalog.likelihoodMethod().name()))
                .isCanonical(catalog.isCanonical())
                .build();
    }

    @NotNull
    public static PurpleQC convert(@NotNull com.hartwig.hmftools.common.purple.PurpleQC purpleQC) {
        return ImmutablePurpleQC.builder()
                .status(ConversionUtil.mapToIterable(purpleQC.status(), PurpleConversion::convert))
                .germlineAberrations(ConversionUtil.mapToIterable(purpleQC.germlineAberrations(), PurpleConversion::convert))
                .amberMeanDepth(purpleQC.amberMeanDepth())
                .contamination(purpleQC.contamination())
                .unsupportedCopyNumberSegments(purpleQC.unsupportedCopyNumberSegments())
                .deletedGenes(purpleQC.deletedGenes())
                .build();
    }

    @NotNull
    public static PurpleAllelicDepth convert(AllelicDepth allelicDepth) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(allelicDepth.alleleReadCount())
                .totalReadCount(allelicDepth.totalReadCount())
                .build();
    }

    @NotNull
    public static PurpleGermlineAberration convert(GermlineAberration i) {
        return PurpleGermlineAberration.valueOf(i.name());
    }

    @NotNull
    public static PurpleQCStatus convert(com.hartwig.hmftools.common.purple.PurpleQCStatus i) {
        return PurpleQCStatus.valueOf(i.name());
    }

    public static PurpleCodingEffect convert(CodingEffect effect) {
        return PurpleCodingEffect.valueOf(effect.name());
    }

    @NotNull
    public static PurpleVariantEffect convert(VariantEffect effect) {
        return PurpleVariantEffect.valueOf(effect.name());
    }
}
