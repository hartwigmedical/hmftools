package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.purple.loader.PurpleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PurpleInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(PurpleInterpreter.class);

    private PurpleInterpreter() {
    }

    @NotNull
    public static PurpleInterpretedData interpret(@NotNull PurpleData purple, @NotNull List<DriverGene> driverGenes,
            @NotNull ChordData chord) {
        List<ReportableVariant> reportableSomaticVariants =
                ReportableVariantFactory.toReportableSomaticVariants(purple.reportableSomaticVariants(), purple.somaticDrivers());
        List<ReportableVariant> additionalSuspectSomaticVariants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(purple.allSomaticVariants(),
                        reportableSomaticVariants,
                        driverGenes);
        LOGGER.info(" Found an additional {} somatic variants that are potentially interesting", additionalSuspectSomaticVariants.size());

        List<ReportableVariant> additionalSuspectGermlineVariants =
                GermlineVariantSelector.selectInterestingUnreportedVariants(purple.allGermlineVariants());
        if (additionalSuspectGermlineVariants != null) {
            LOGGER.info(" Found an additional {} germline variants that are potentially interesting",
                    additionalSuspectGermlineVariants.size());
        }

        List<GainLoss> nearReportableSomaticGains = CopyNumberSelector.selectNearReportableSomaticGains(purple.allSomaticGeneCopyNumbers(),
                purple.purityContext().bestFit().ploidy(),
                purple.reportableSomaticGainsLosses(),
                driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting",
                nearReportableSomaticGains.size());

        List<GainLoss> additionalSuspectSomaticGainsLosses =
                CopyNumberSelector.selectInterestingUnreportedGainsLosses(purple.allSomaticGainsLosses(),
                        purple.reportableSomaticGainsLosses());
        LOGGER.info(" Found an additional {} somatic gains/losses that are potentially interesting",
                additionalSuspectSomaticGainsLosses.size());

        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(purple.allSomaticGeneCopyNumbers(),
                        purple.purityContext().microsatelliteStatus(),
                        chord.hrStatus());
        LOGGER.info(" Found an additional {} suspect gene copy numbers with LOH", suspectGeneCopyNumbersWithLOH.size());

        List<ReportableVariant> reportableGermlineVariants = purple.reportableGermlineVariants() != null
                ? ReportableVariantFactory.toReportableGermlineVariants(purple.reportableGermlineVariants(), purple.germlineDrivers())
                : null;

        return ImmutablePurpleInterpretedData.builder()
                .fit(createFit(purple))
                .characteristics(createCharacteristics(purple))
                .allSomaticVariants(purple.allSomaticVariants())
                .reportableSomaticVariants(reportableSomaticVariants)
                .additionalSuspectSomaticVariants(additionalSuspectSomaticVariants)
                .allGermlineVariants(purple.allGermlineVariants())
                .reportableGermlineVariants(reportableGermlineVariants)
                .additionalSuspectGermlineVariants(additionalSuspectGermlineVariants)
                .allSomaticGeneCopyNumbers(purple.allSomaticGeneCopyNumbers())
                .suspectGeneCopyNumbersWithLOH(suspectGeneCopyNumbersWithLOH)
                .allSomaticGainsLosses(purple.allSomaticGainsLosses())
                .reportableSomaticGainsLosses(purple.reportableSomaticGainsLosses())
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsLosses(additionalSuspectSomaticGainsLosses)
                .allGermlineDeletions(purple.allGermlineDeletions())
                .reportableGermlineDeletions(purple.reportableGermlineDeletions())
                .build();
    }

    @NotNull
    private static PurityPloidyFit createFit(@NotNull PurpleData purple) {
        return ImmutablePurityPloidyFit.builder()
                .qc(purple.purityContext().qc())
                .hasReliableQuality(purple.purityContext().qc().pass())
                .fittedPurityMethod(purple.purityContext().method())
                .hasReliablePurity(PurityContext.checkHasReliablePurity(purple.purityContext()))
                .purity(purple.purityContext().bestFit().purity())
                .minPurity(purple.purityContext().score().minPurity())
                .maxPurity(purple.purityContext().score().maxPurity())
                .ploidy(purple.purityContext().bestFit().ploidy())
                .minPloidy(purple.purityContext().score().minPloidy())
                .maxPloidy(purple.purityContext().score().maxPloidy())
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createCharacteristics(@NotNull PurpleData purple) {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(purple.purityContext().wholeGenomeDuplication())
                .microsatelliteIndelsPerMb(purple.purityContext().microsatelliteIndelsPerMb())
                .microsatelliteStatus(purple.purityContext().microsatelliteStatus())
                .tumorMutationalBurdenPerMb(purple.purityContext().tumorMutationalBurdenPerMb())
                .tumorMutationalLoad(purple.purityContext().tumorMutationalLoad())
                .tumorMutationalLoadStatus(purple.purityContext().tumorMutationalLoadStatus())
                .svTumorMutationalBurden(purple.purityContext().svTumorMutationalBurden())
                .build();
    }
}
