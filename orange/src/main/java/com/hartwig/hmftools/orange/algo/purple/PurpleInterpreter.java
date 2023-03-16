package com.hartwig.hmftools.orange.algo.purple;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.*;
import com.hartwig.hmftools.common.drivercatalog.panel.*;
import com.hartwig.hmftools.common.purple.*;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.purple.*;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import com.hartwig.hmftools.orange.algo.linx.BreakendUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

public class PurpleInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(PurpleInterpreter.class);

    private static final int MAX_LENGTH_FOR_IMPLIED_DELS = 1500;

    @NotNull
    private final PurpleVariantFactory purpleVariantFactory;
    @NotNull
    private final GermlineGainLossFactory germlineGainLossFactory;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final LinxRecord linx;
    @Nullable
    private final ChordRecord chord;

    public PurpleInterpreter(@NotNull final PurpleVariantFactory purpleVariantFactory,
                             @NotNull final GermlineGainLossFactory germlineGainLossFactory, @NotNull final List<DriverGene> driverGenes,
                             @NotNull final LinxRecord linx, @Nullable final ChordRecord chord) {
        this.purpleVariantFactory = purpleVariantFactory;
        this.germlineGainLossFactory = germlineGainLossFactory;
        this.driverGenes = driverGenes;
        this.linx = linx;
        this.chord = chord;
    }

    @NotNull
    public PurpleRecord interpret(@NotNull PurpleData purple) {
        LOGGER.info("Analysing purple data");
        List<PurpleVariant> allSomaticVariants = purpleVariantFactory.create(purple.allSomaticVariants());
        List<PurpleVariant> reportableSomaticVariants = purpleVariantFactory.create(purple.reportableSomaticVariants());
        List<PurpleVariant> additionalSuspectSomaticVariants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(allSomaticVariants, reportableSomaticVariants, driverGenes);
        LOGGER.info(" Found an additional {} somatic variants that are potentially interesting", additionalSuspectSomaticVariants.size());

        List<PurpleVariant> allGermlineVariants = purpleVariantFactory.create(purple.allGermlineVariants());
        List<PurpleVariant> reportableGermlineVariants = purpleVariantFactory.create(purple.reportableGermlineVariants());
        List<PurpleVariant> additionalSuspectGermlineVariants =
                GermlineVariantSelector.selectInterestingUnreportedVariants(allGermlineVariants);
        if (additionalSuspectGermlineVariants != null) {
            LOGGER.info(" Found an additional {} germline variants that are potentially interesting",
                    additionalSuspectGermlineVariants.size());
        }

        List<PurpleGainLoss> allSomaticGainsLosses = extractAllGainsLosses(purple.purityContext().qc().status(),
                purple.purityContext().bestFit().ploidy(),
                purple.purityContext().targeted(),
                purple.allSomaticGeneCopyNumbers());
        List<PurpleGainLoss> reportableSomaticGainsLosses = somaticGainsLossesFromDrivers(purple.somaticDrivers());

        List<PurpleGainLoss> nearReportableSomaticGains =
                CopyNumberSelector.selectNearReportableSomaticGains(purple.allSomaticGeneCopyNumbers(),
                        purple.purityContext().bestFit().ploidy(),
                        allSomaticGainsLosses,
                        driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting",
                nearReportableSomaticGains.size());

        List<PurpleGainLoss> additionalSuspectSomaticGainsLosses =
                CopyNumberSelector.selectInterestingUnreportedGainsLosses(allSomaticGainsLosses, reportableSomaticGainsLosses);
        LOGGER.info(" Found an additional {} somatic gains/losses that are potentially interesting",
                additionalSuspectSomaticGainsLosses.size());

        List<GermlineDeletion> allGermlineDeletions = purple.allGermlineDeletions();
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(purple.allSomaticGeneCopyNumbers(),
                        allGermlineDeletions,
                        purple.purityContext().microsatelliteStatus(),
                        chord != null ? chord.hrStatus() : null);
        LOGGER.info(" Found an additional {} suspect gene copy numbers with LOH", suspectGeneCopyNumbersWithLOH.size());

        List<PurpleGainLoss> allGermlineFullLosses = null;
        List<PurpleGainLoss> reportableGermlineFullLosses = null;

        if (allGermlineDeletions != null) {
            List<GermlineDeletion> impliedDeletions = implyDeletionsFromBreakends(allGermlineDeletions,
                    linx.reportableGermlineBreakends(),
                    purple.allGermlineStructuralVariants(),
                    linx.allGermlineStructuralVariants());
            LOGGER.info(" Implied {} additional reportable germline deletions from breakends", impliedDeletions.size());

            List<GermlineDeletion> mergedGermlineDeletions = Lists.newArrayList();
            mergedGermlineDeletions.addAll(allGermlineDeletions);
            mergedGermlineDeletions.addAll(impliedDeletions);

            Map<PurpleGainLoss, GermlineDeletion> deletionMap =
                    germlineGainLossFactory.mapDeletions(mergedGermlineDeletions, purple.allSomaticGeneCopyNumbers());

            allGermlineFullLosses = Lists.newArrayList(deletionMap.keySet());
            reportableGermlineFullLosses = selectReportable(deletionMap);

            LOGGER.info(" Resolved {} germline losses of which {} are reportable",
                    allGermlineFullLosses.size(),
                    reportableGermlineFullLosses.size());
        }

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purple))
                .characteristics(createCharacteristics(purple))
                .somaticDrivers(() -> purple.somaticDrivers().stream().map(PurpleInterpreter::asPurpleDriver).iterator())
                .germlineDrivers(() -> Objects.requireNonNullElseGet(purple.germlineDrivers(), List::<DriverCatalog>of).stream().map(PurpleInterpreter::asPurpleDriver).iterator())
                .allSomaticVariants(allSomaticVariants)
                .reportableSomaticVariants(reportableSomaticVariants)
                .additionalSuspectSomaticVariants(additionalSuspectSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .additionalSuspectGermlineVariants(additionalSuspectGermlineVariants)
                .allSomaticCopyNumbers(() -> purple.allSomaticCopyNumbers().stream().map(PurpleInterpreter::asPurpleCopyNumber).iterator())
                .allSomaticGeneCopyNumbers(() -> purple.allSomaticGeneCopyNumbers().stream().map(PurpleInterpreter::asPurpleGeneCopyNumber).iterator())
                .suspectGeneCopyNumbersWithLOH(() -> suspectGeneCopyNumbersWithLOH.stream().map(PurpleInterpreter::asPurpleGeneCopyNumber).iterator())
                .allSomaticGainsLosses(allSomaticGainsLosses)
                .reportableSomaticGainsLosses(reportableSomaticGainsLosses)
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsLosses(additionalSuspectSomaticGainsLosses)
                .reportableGermlineFullLosses(reportableGermlineFullLosses)
                .build();
    }

    private static PurpleCopyNumber asPurpleCopyNumber(com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .averageTumorCopyNumber(copyNumber.averageTumorCopyNumber())
                .build();
    }

    private static PurpleGeneCopyNumber asPurpleGeneCopyNumber(GeneCopyNumber geneCopyNumber) {
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

    @NotNull
    @VisibleForTesting
    static List<GermlineDeletion> implyDeletionsFromBreakends(@NotNull List<GermlineDeletion> allGermlineDeletions,
            @Nullable List<LinxBreakend> reportableGermlineBreakends, @NotNull List<StructuralVariant> allPurpleGermlineSvs,
            @Nullable List<LinxSvAnnotation> allLinxGermlineSvAnnotations) {
        if (reportableGermlineBreakends == null || allLinxGermlineSvAnnotations == null) {
            LOGGER.warn("Linx germline data is missing while purple germline data is present!");
            return Lists.newArrayList();
        }

        List<GermlineDeletion> impliedDeletions = Lists.newArrayList();
        for (Pair<LinxBreakend, LinxBreakend> breakendPair : BreakendUtil.createPairsPerSvId(reportableGermlineBreakends)) {
            LinxBreakend first = breakendPair.getLeft();
            LinxBreakend second = breakendPair.getRight();

            boolean bothReported = first.reportedDisruption() && second.reportedDisruption();
            boolean bothDel = first.type() == LinxBreakendType.DEL && second.type() == LinxBreakendType.DEL;
            boolean sameGene = first.gene().equals(second.gene());
            boolean sameTranscript = first.transcriptId().equals(second.transcriptId());
            boolean noWildTypeRemaining = first.undisruptedCopyNumber() < 0.5 && second.undisruptedCopyNumber() < 0.5;

            StructuralVariant sv = findBySvId(allPurpleGermlineSvs, allLinxGermlineSvAnnotations, first.svId());
            boolean meetsMaxLength = false;
            if (sv != null) {
                meetsMaxLength = Math.abs(sv.start().position() - sv.end().position()) <= MAX_LENGTH_FOR_IMPLIED_DELS;
            }

            boolean hasNoExistingGermlineDel = !hasGermlineDeletionInGene(allGermlineDeletions, first.gene());

            if (bothReported && bothDel && sameGene && sameTranscript && noWildTypeRemaining && meetsMaxLength
                    && hasNoExistingGermlineDel) {
                impliedDeletions.add(new GermlineDeletion(first.gene(),
                        first.chromosome(),
                        first.chrBand(),
                        0,
                        0,
                        0,
                        0,
                        0,
                        GermlineDetectionMethod.SEGMENT,
                        GermlineStatus.HET_DELETION,
                        GermlineStatus.HOM_DELETION,
                        1D,
                        0D,
                        Strings.EMPTY,
                        0,
                        true));
            }
        }
        return impliedDeletions;
    }

    @Nullable
    private static StructuralVariant findBySvId(@NotNull List<StructuralVariant> allPurpleSvs,
            @NotNull List<LinxSvAnnotation> allLinxSvAnnotations, int svIdToFind) {
        LinxSvAnnotation match = null;
        int index = 0;
        while (match == null && index < allLinxSvAnnotations.size()) {
            LinxSvAnnotation svAnnotation = allLinxSvAnnotations.get(index);
            if (svAnnotation.svId() == svIdToFind) {
                match = svAnnotation;
            }
            index++;
        }

        if (match != null) {
            for (StructuralVariant structuralVariant : allPurpleSvs) {
                if (structuralVariant.id().equals(match.vcfId())) {
                    return structuralVariant;
                }
            }
        }

        LOGGER.warn("Could not find germline structural variant with svId: {}", svIdToFind);
        return null;
    }

    private static boolean hasGermlineDeletionInGene(@NotNull List<GermlineDeletion> germlineDeletions, @NotNull String geneToFind) {
        for (GermlineDeletion deletion : germlineDeletions) {
            if (deletion.GeneName.equals(geneToFind)) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<PurpleGainLoss> selectReportable(@NotNull Map<PurpleGainLoss, GermlineDeletion> deletionMap) {
        List<PurpleGainLoss> reportable = Lists.newArrayList();
        for (Map.Entry<PurpleGainLoss, GermlineDeletion> entry : deletionMap.entrySet()) {
            PurpleGainLoss gainLoss = entry.getKey();
            GermlineDeletion deletion = entry.getValue();
            if (deletion.Reported) {
                reportable.add(gainLoss);
            }
        }
        return reportable;
    }

    @NotNull
    private static List<PurpleGainLoss> extractAllGainsLosses(@NotNull Set<PurpleQCStatus> qcStatus, double ploidy, boolean isTargetRegions,
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
                    .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                    .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
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
    private static List<PurpleGainLoss> somaticGainsLossesFromDrivers(@NotNull List<DriverCatalog> drivers) {
        List<PurpleGainLoss> gainsLosses = Lists.newArrayList();

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
    private static PurpleGainLoss toGainLoss(@NotNull DriverCatalog driver) {
        return ImmutablePurpleGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(CopyNumberInterpretationUtil.fromCNADriver(driver))
                .minCopies(Math.round(Math.max(0, driver.minCopyNumber())))
                .maxCopies(Math.round(Math.max(0, driver.maxCopyNumber())))
                .build();
    }

    @NotNull
    private static PurpleFit createFit(@NotNull PurpleData purple) {
        return ImmutablePurpleFit.builder()
                .qc(createQC(purple.purityContext().qc()))
                .hasSufficientQuality(purple.purityContext().qc().pass())
                .fittedPurityMethod(PurpleFittedPurityMethod.valueOf(purple.purityContext().method().name()))
                .containsTumorCells(purple.purityContext().method() != FittedPurityMethod.NO_TUMOR)
                .purity(purple.purityContext().bestFit().purity())
                .minPurity(purple.purityContext().score().minPurity())
                .maxPurity(purple.purityContext().score().maxPurity())
                .ploidy(purple.purityContext().bestFit().ploidy())
                .minPloidy(purple.purityContext().score().minPloidy())
                .maxPloidy(purple.purityContext().score().maxPloidy())
                .build();
    }

    private static PurpleQC createQC(@NotNull com.hartwig.hmftools.common.purple.PurpleQC purpleQC) {
        return ImmutablePurpleQC.builder()
                .status(() -> purpleQC.status().stream().map(i -> com.hartwig.hmftools.datamodel.purple.PurpleQCStatus.valueOf(i.name())).iterator())
                .germlineAberrations(() -> purpleQC.germlineAberrations().stream().map(i -> PurpleGermlineAberration.valueOf(i.name())).iterator())
                .amberMeanDepth(purpleQC.amberMeanDepth())
                .contamination(purpleQC.contamination())
                .unsupportedCopyNumberSegments(purpleQC.unsupportedCopyNumberSegments())
                .deletedGenes(purpleQC.deletedGenes())
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createCharacteristics(@NotNull PurpleData purple) {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(purple.purityContext().wholeGenomeDuplication())
                .microsatelliteIndelsPerMb(purple.purityContext().microsatelliteIndelsPerMb())
                .microsatelliteStatus(PurpleMicrosatelliteStatus.valueOf(purple.purityContext().microsatelliteStatus().name()))
                .tumorMutationalBurdenPerMb(purple.purityContext().tumorMutationalBurdenPerMb())
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext().tumorMutationalBurdenStatus().name()))
                .tumorMutationalLoad(purple.purityContext().tumorMutationalLoad())
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext().tumorMutationalLoadStatus().name()))
                .svTumorMutationalBurden(purple.purityContext().svTumorMutationalBurden())
                .build();
    }
}
