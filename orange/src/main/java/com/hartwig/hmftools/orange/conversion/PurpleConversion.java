package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.*;
import org.jetbrains.annotations.NotNull;

public class PurpleConversion {

    private PurpleConversion() {
    }

    public static PurpleCopyNumber asPurpleCopyNumber(com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .averageTumorCopyNumber(copyNumber.averageTumorCopyNumber())
                .build();
    }

    public static PurpleGeneCopyNumber asPurpleGeneCopyNumber(GeneCopyNumber geneCopyNumber) {
        return ImmutablePurpleGeneCopyNumber.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .geneName(geneCopyNumber.geneName())
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .minMinorAlleleCopyNumber(geneCopyNumber.minMinorAlleleCopyNumber())
                .build();
    }

    public static PurpleDriver asPurpleDriver(DriverCatalog catalog) {
        return ImmutablePurpleDriver.builder()
                .gene(catalog.gene())
                .transcript(catalog.transcript())
                .driver(PurpleDriverType.valueOf(catalog.driver().name()))
                .driverLikelihood(catalog.driverLikelihood())
                .likelihoodMethod(PurpleLikelihoodMethod.valueOf(catalog.likelihoodMethod().name()))
                .isCanonical(catalog.isCanonical())
                .build();
    }

    public static PurpleQC createQC(@NotNull com.hartwig.hmftools.common.purple.PurpleQC purpleQC) {
        return ImmutablePurpleQC.builder()
                .status(() -> purpleQC.status().stream().map(i -> PurpleQCStatus.valueOf(i.name())).iterator())
                .germlineAberrations(() -> purpleQC.germlineAberrations().stream().map(i -> PurpleGermlineAberration.valueOf(i.name())).iterator())
                .amberMeanDepth(purpleQC.amberMeanDepth())
                .contamination(purpleQC.contamination())
                .unsupportedCopyNumberSegments(purpleQC.unsupportedCopyNumberSegments())
                .deletedGenes(purpleQC.deletedGenes())
                .build();
    }
}
