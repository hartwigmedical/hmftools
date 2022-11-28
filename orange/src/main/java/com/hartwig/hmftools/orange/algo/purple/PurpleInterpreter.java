package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.drivercatalog.AmplificationDrivers;
import com.hartwig.hmftools.common.drivercatalog.DeletionDrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMap;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGenePanel;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;

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
        List<GainLoss> allSomaticGainsLosses = extractAllGainsLosses(purple.purityContext().qc().status(),
                purple.purityContext().bestFit().ploidy(),
                purple.purityContext().targeted(),
                purple.allSomaticGeneCopyNumbers());
        List<GainLoss> reportableSomaticGainsLosses = somaticGainsLossesFromDrivers(purple.somaticDrivers());

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
                allSomaticGainsLosses,
                driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting",
                nearReportableSomaticGains.size());

        List<GainLoss> additionalSuspectSomaticGainsLosses =
                CopyNumberSelector.selectInterestingUnreportedGainsLosses(allSomaticGainsLosses,
                        reportableSomaticGainsLosses);
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
                .allSomaticGainsLosses(allSomaticGainsLosses)
                .reportableSomaticGainsLosses(reportableSomaticGainsLosses)
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsLosses(additionalSuspectSomaticGainsLosses)
                .allGermlineDeletions(purple.allGermlineDeletions())
                .reportableGermlineDeletions(purple.reportableGermlineDeletions())
                .build();
    }

    @NotNull
    private static List<GainLoss> extractAllGainsLosses(@NotNull Set<PurpleQCStatus> qcStatus, double ploidy, boolean isTargetRegions,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers) {
        List<DriverGene> allGenes = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) {
            allGenes.add(ImmutableDriverGene.builder()
                    .gene(geneCopyNumber.geneName())
                    .reportMissenseAndInframe(false)
                    .reportNonsenseAndFrameshift(false)
                    .reportSplice(false)
                    .reportDeletion(true)
                    .reportDisruption(false)
                    .reportAmplification(true)
                    .reportSomaticHotspot(false)
                    .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                    .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                    .reportGermlineDisruption(false)
                    .likelihoodType(DriverCategory.ONCO)
                    .reportPGX(false)
                    .build());
        }

        DriverGenePanel allGenesPanel = ImmutableDriverGenePanel.builder().driverGenes(allGenes).build();
        AmplificationDrivers ampDrivers = new AmplificationDrivers(qcStatus, allGenesPanel);
        DeletionDrivers delDrivers = new DeletionDrivers(qcStatus, allGenesPanel);

        List<DriverCatalog> allGainLosses = Lists.newArrayList();
        allGainLosses.addAll(ampDrivers.amplifications(ploidy, allGeneCopyNumbers, isTargetRegions));
        allGainLosses.addAll(delDrivers.deletions(allGeneCopyNumbers, isTargetRegions));

        return somaticGainsLossesFromDrivers(allGainLosses);
    }

    @NotNull
    private static List<GainLoss> somaticGainsLossesFromDrivers(@NotNull List<DriverCatalog> drivers) {
        List<GainLoss> gainsLosses = Lists.newArrayList();

        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(drivers);
        for (DriverCatalogKey key : geneDriverMap.keySet()) {
            DriverCatalog geneDriver = geneDriverMap.get(key);

            if (geneDriver.driver() == DriverType.AMP || geneDriver.driver() == DriverType.PARTIAL_AMP
                    || geneDriver.driver() == DriverType.DEL) {
                gainsLosses.add(toGainLoss(geneDriver));
            }
        }
        return gainsLosses;
    }

    @NotNull
    private static GainLoss toGainLoss(@NotNull DriverCatalog driver) {
        return ImmutableGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                .minCopies(Math.round(Math.max(0, driver.minCopyNumber())))
                .maxCopies(Math.round(Math.max(0, driver.maxCopyNumber())))
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
