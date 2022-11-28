package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PurpleTestFactory {

    private PurpleTestFactory() {
    }

    @NotNull
    public static PurpleData createMinimalTestPurpleData()
    {
        PurityContext minimalContext = ImmutablePurityContext.builder()
                .version(Strings.EMPTY)
                .gender(Gender.FEMALE)
                .runMode(RunMode.TUMOR_GERMLINE)
                .targeted(false)
                .bestFit(emptyFit())
                .method(FittedPurityMethod.NORMAL)
                .score(emptyScore())
                .qc(qcPass())
                .polyClonalProportion(0D)
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(0D)
                .tumorMutationalBurdenPerMb(0D)
                .tumorMutationalLoad(0)
                .svTumorMutationalBurden(0)
                .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.UNKNOWN)
                .build();

        return ImmutablePurpleData.builder().purityContext(minimalContext).build();
    }

    @NotNull
    private static FittedPurity emptyFit()
    {
        return ImmutableFittedPurity.builder()
                .purity(0D)
                .normFactor(0D)
                .ploidy(0D)
                .score(0D)
                .diploidProportion(0D)
                .somaticPenalty(0D)
                .build();
    }

    @NotNull
    private static FittedPurityScore emptyScore()
    {
        return ImmutableFittedPurityScore.builder()
                .minPurity(0D)
                .maxPurity(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .minDiploidProportion(0D)
                .maxDiploidProportion(0D)
                .build();
    }

    @NotNull
    private static PurpleQC qcPass()
    {
        return ImmutablePurpleQC.builder()
                .method(FittedPurityMethod.NORMAL)
                .amberMeanDepth(0)
                .copyNumberSegments(1)
                .unsupportedCopyNumberSegments(0)
                .deletedGenes(0)
                .purity(0.5)
                .contamination(0D)
                .cobaltGender(Gender.FEMALE)
                .amberGender(Gender.FEMALE)
                .build();
    }
}
