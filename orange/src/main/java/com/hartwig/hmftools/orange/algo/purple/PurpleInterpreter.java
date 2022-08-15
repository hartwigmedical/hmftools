package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.loader.PurpleData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PurpleInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(PurpleInterpreter.class);

    private PurpleInterpreter() {
    }

    @NotNull
    public static PurpleInterpretedData interpret(@NotNull PurpleData purple, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<DriverGene> driverGenes, @NotNull ChordData chord) {
        List<ReportableVariant> additionalSuspectSomaticVariants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(purple.allSomaticVariants(),
                        purple.reportableSomaticVariants(),
                        evidences,
                        driverGenes);
        LOGGER.info(" Found an additional {} somatic variants that are potentially interesting", additionalSuspectSomaticVariants.size());

        List<ReportableVariant> additionalSuspectGermlineVariants =
                GermlineVariantSelector.selectInterestingUnreportedVariants(purple.allGermlineVariants());
        LOGGER.info(" Found an additional {} germline variants that are potentially interesting", additionalSuspectGermlineVariants.size());

        List<GainLoss> nearReportableSomaticGains = CopyNumberSelector.selectNearReportableSomaticGains(purple.allSomaticGeneCopyNumbers(),
                purple.ploidy(),
                purple.reportableSomaticGainsLosses(),
                driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting",
                nearReportableSomaticGains.size());

        List<GainLoss> additionalSuspectSomaticGainsLosses =
                CopyNumberSelector.selectInterestingUnreportedGainsLosses(purple.allSomaticGainsLosses(),
                        purple.reportableSomaticGainsLosses(),
                        evidences);
        LOGGER.info(" Found an additional {} somatic gains/losses that are potentially interesting",
                additionalSuspectSomaticGainsLosses.size());

        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(purple.allSomaticGeneCopyNumbers(),
                        purple.microsatelliteStatus(),
                        chord.hrStatus());
        LOGGER.info(" Found an additional {} suspect gene copy numbers with LOH", suspectGeneCopyNumbersWithLOH.size());

        return ImmutablePurpleInterpretedData.builder()
                .fit(createFit(purple))
                .characteristics(createCharacteristics(purple))
                .allSomaticVariants(purple.allSomaticVariants())
                .reportableSomaticVariants(purple.reportableSomaticVariants())
                .additionalSuspectSomaticVariants(additionalSuspectSomaticVariants)
                .allGermlineVariants(purple.allGermlineVariants())
                .reportableGermlineVariants(purple.reportableGermlineVariants())
                .additionalSuspectGermlineVariants(additionalSuspectGermlineVariants)
                .allSomaticGeneCopyNumbers(purple.allSomaticGeneCopyNumbers())
                .suspectGeneCopyNumbersWithLOH(suspectGeneCopyNumbersWithLOH)
                .allSomaticGainsLosses(purple.allSomaticGainsLosses())
                .reportableSomaticGainsLosses(purple.reportableSomaticGainsLosses())
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsLosses(additionalSuspectSomaticGainsLosses)
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
