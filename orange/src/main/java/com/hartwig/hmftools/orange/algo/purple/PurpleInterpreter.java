package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.purple.PurpleData;

import org.jetbrains.annotations.NotNull;

public final class PurpleInterpreter {

    private PurpleInterpreter() {
    }

    @NotNull
    public static PurpleInterpretedData interpret(@NotNull PurpleData purple) {
        return ImmutablePurpleInterpretedData.builder()
                .fit(createFit(purple))
                .characteristics(createCharacteristics(purple))
                .allSomaticVariants(purple.allSomaticVariants())
                .reportableSomaticVariants(purple.reportableSomaticVariants())
                .allGermlineVariants(purple.allGermlineVariants())
                .reportableGermlineVariants(purple.reportableGermlineVariants())
                .allSomaticGeneCopyNumbers(purple.allSomaticGeneCopyNumbers())
                .allSomaticGainsLosses(purple.allSomaticGainsLosses())
                .reportableSomaticGainsLosses(purple.reportableSomaticGainsLosses())
                .allGermlineDeletions(purple.allGermlineDeletions())
                .reportableGermlineDeletions(purple.reportableGermlineDeletions())
                .copyNumberPerChromosome(purple.copyNumberPerChromosome())
                .build();
    }

    @NotNull
    private static PurityPloidyFit createFit(@NotNull PurpleData purple) {
        return ImmutablePurityPloidyFit.builder()
                .qc(purple.qc())
                .hasReliableQuality(purple.hasReliableQuality())
                .fittedPurityMethod(purple.fittedPurityMethod())
                .hasReliablePurity(purple.hasReliablePurity())
                .purity(purple.purity())
                .minPurity(purple.minPurity())
                .maxPurity(purple.maxPurity())
                .ploidy(purple.ploidy())
                .minPloidy(purple.minPloidy())
                .maxPloidy(purple.maxPloidy())
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createCharacteristics(@NotNull PurpleData purple) {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(purple.wholeGenomeDuplication())
                .microsatelliteIndelsPerMb(purple.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purple.microsatelliteStatus())
                .tumorMutationalBurdenPerMb(purple.tumorMutationalBurdenPerMb())
                .tumorMutationalLoad(purple.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purple.tumorMutationalLoadStatus())
                .svTumorMutationalBurden(purple.svTumorMutationalBurden())
                .build();
    }
}
