package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
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
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.linx.BreakendUtil;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    private static final int MAX_LENGTH_FOR_IMPLIED_DELS = 1500;
    @NotNull
    private final PurpleVariantFactory purpleVariantFactory;
    @NotNull
    private final GermlineGainLossFactory germlineGainLossFactory;
    @NotNull
    private final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final LinxRecord linx;
    @Nullable
    private final ChordData chord;

    public PurpleInterpreter(@NotNull final PurpleVariantFactory purpleVariantFactory,
            @NotNull final GermlineGainLossFactory germlineGainLossFactory,
            @NotNull final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory,
            @NotNull final List<DriverGene> driverGenes, @NotNull final LinxRecord linx, @Nullable final ChordData chord)
    {
        this.purpleVariantFactory = purpleVariantFactory;
        this.germlineGainLossFactory = germlineGainLossFactory;
        this.germlineLossOfHeterozygosityFactory = germlineLossOfHeterozygosityFactory;
        this.driverGenes = driverGenes;
        this.linx = linx;
        this.chord = chord;
    }

    @NotNull
    public PurpleRecord interpret(@NotNull PurpleData purple)
    {
        LOGGER.info("Analysing purple data");

        List<PurpleVariant> allSomaticVariants = purpleVariantFactory.fromPurpleVariantContext(purple.allSomaticVariants());

        List<PurpleVariant> reportableSomaticVariants = purpleVariantFactory.fromPurpleVariantContext(purple.reportableSomaticVariants());
        List<PurpleVariant> additionalSuspectSomaticVariants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(allSomaticVariants, reportableSomaticVariants, driverGenes);
        LOGGER.info(" Found an additional {} somatic variants that are potentially interesting", additionalSuspectSomaticVariants.size());

        List<PurpleVariant> allGermlineVariants = purpleVariantFactory.fromPurpleVariantContext(purple.allGermlineVariants());
        List<PurpleVariant> reportableGermlineVariants = purpleVariantFactory.fromPurpleVariantContext(purple.reportableGermlineVariants());
        List<PurpleVariant> additionalSuspectGermlineVariants =
                GermlineVariantSelector.selectInterestingUnreportedVariants(allGermlineVariants);
        if(additionalSuspectGermlineVariants != null)
        {
            LOGGER.info(" Found an additional {} germline variants that are potentially interesting", additionalSuspectGermlineVariants.size());
        }

        List<PurpleGainLoss> allSomaticGainsLosses =
                extractAllGainsLosses(purple.purityContext().qc().status(), purple.purityContext().gender(), purple.purityContext()
                        .bestFit()
                        .ploidy(), purple.purityContext().targeted(), purple.allSomaticGeneCopyNumbers());
        List<PurpleGainLoss> reportableSomaticGainsLosses = somaticGainsLossesFromDrivers(purple.somaticDrivers());

        List<PurpleGainLoss> nearReportableSomaticGains =
                CopyNumberSelector.selectNearReportableSomaticGains(purple.allSomaticGeneCopyNumbers(), purple.purityContext()
                        .bestFit()
                        .ploidy(), allSomaticGainsLosses, driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting", nearReportableSomaticGains.size());

        List<PurpleGainLoss> additionalSuspectSomaticGainsLosses =
                CopyNumberSelector.selectInterestingUnreportedGainsLosses(allSomaticGainsLosses, reportableSomaticGainsLosses);
        LOGGER.info(" Found an additional {} somatic gains/losses that are potentially interesting", additionalSuspectSomaticGainsLosses.size());

        List<GermlineDeletion> allGermlineDeletions = purple.allGermlineDeletions();
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(purple.allSomaticGeneCopyNumbers(), allGermlineDeletions, purple.purityContext()
                        .microsatelliteStatus(), chord != null ? chord.hrStatus() : null);
        LOGGER.info(" Found an additional {} suspect gene copy numbers with LOH", suspectGeneCopyNumbersWithLOH.size());

        List<PurpleGainLoss> allGermlineFullLosses = null;
        List<PurpleGainLoss> reportableGermlineFullLosses = null;
        List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities = null;
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = null;

        if(allGermlineDeletions != null)
        {
            List<GermlineDeletion> impliedDeletions =
                    implyDeletionsFromBreakends(allGermlineDeletions, linx.reportableGermlineBreakends(), purple.allGermlineStructuralVariants(), linx.allGermlineStructuralVariants(), driverGenes);
            LOGGER.info(" Implied {} additional reportable germline deletions from breakends", impliedDeletions.size());

            List<GermlineDeletion> mergedGermlineDeletions = Lists.newArrayList();
            mergedGermlineDeletions.addAll(allGermlineDeletions);
            mergedGermlineDeletions.addAll(impliedDeletions);

            Map<PurpleGainLoss, Boolean> fullLossToReportability =
                    germlineGainLossFactory.getReportabilityMap(mergedGermlineDeletions, purple.allSomaticGeneCopyNumbers());

            allGermlineFullLosses = Lists.newArrayList(fullLossToReportability.keySet());
            reportableGermlineFullLosses = selectReportableGainLosses(fullLossToReportability);

            LOGGER.info(" Resolved {} germline losses of which {} are reportable", allGermlineFullLosses.size(), reportableGermlineFullLosses.size());

            Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability =
                    germlineLossOfHeterozygosityFactory.getReportabilityMap(mergedGermlineDeletions, purple.allSomaticGeneCopyNumbers());

            allGermlineLossOfHeterozygosities = Lists.newArrayList(lossOfHeterozygosityToReportability.keySet());
            reportableGermlineLossOfHeterozygosities = selectReportableLossOfHeterozygosities(lossOfHeterozygosityToReportability);

            LOGGER.info(" Resolved {} germline heterozygous deletions of which {} are reportable", allGermlineLossOfHeterozygosities.size(), reportableGermlineLossOfHeterozygosities.size());
        }

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purple))
                .characteristics(createCharacteristics(purple))
                .somaticDrivers(ConversionUtil.mapToIterable(purple.somaticDrivers(), PurpleConversion::convert))
                .germlineDrivers(ConversionUtil.mapToIterable(purple.germlineDrivers(), PurpleConversion::convert))
                .allSomaticVariants(allSomaticVariants)
                .reportableSomaticVariants(reportableSomaticVariants)
                .additionalSuspectSomaticVariants(additionalSuspectSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .additionalSuspectGermlineVariants(additionalSuspectGermlineVariants)
                .allSomaticCopyNumbers(ConversionUtil.mapToIterable(purple.allSomaticCopyNumbers(), PurpleConversion::convert))
                .allSomaticGeneCopyNumbers(ConversionUtil.mapToIterable(purple.allSomaticGeneCopyNumbers(), PurpleConversion::convert))
                .suspectGeneCopyNumbersWithLOH(ConversionUtil.mapToIterable(suspectGeneCopyNumbersWithLOH, PurpleConversion::convert))
                .allSomaticGainsLosses(allSomaticGainsLosses)
                .reportableSomaticGainsLosses(reportableSomaticGainsLosses)
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsLosses(additionalSuspectSomaticGainsLosses)
                .allGermlineDeletions(ConversionUtil.mapToIterable(purple.allGermlineDeletions(), PurpleConversion::convert))
                .allGermlineFullLosses(allGermlineFullLosses)
                .reportableGermlineFullLosses(reportableGermlineFullLosses)
                .allGermlineLossOfHeterozygosities(allGermlineLossOfHeterozygosities)
                .reportableGermlineLossOfHeterozygosities(reportableGermlineLossOfHeterozygosities)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<GermlineDeletion> implyDeletionsFromBreakends(@NotNull List<GermlineDeletion> allGermlineDeletions,
            @Nullable List<LinxBreakend> reportableGermlineBreakends, @NotNull List<StructuralVariant> allPurpleGermlineSvs,
            @Nullable List<LinxSvAnnotation> allLinxGermlineSvAnnotations, @NotNull List<DriverGene> driverGenes)
    {
        if(reportableGermlineBreakends == null || allLinxGermlineSvAnnotations == null)
        {
            LOGGER.warn("Linx germline data is missing while purple germline data is present!");
            return Lists.newArrayList();
        }

        List<GermlineDeletion> impliedDeletions = Lists.newArrayList();
        for(Pair<LinxBreakend, LinxBreakend> breakendPair : BreakendUtil.createPairsPerSvId(reportableGermlineBreakends))
        {
            LinxBreakend first = breakendPair.getLeft();
            LinxBreakend second = breakendPair.getRight();

            boolean bothReported = first.reported() && second.reported();
            boolean bothDel = first.type() == LinxBreakendType.DEL && second.type() == LinxBreakendType.DEL;
            boolean sameGene = first.gene().equals(second.gene());
            boolean sameTranscript = first.transcript().equals(second.transcript());

            StructuralVariant sv = findBySvId(allPurpleGermlineSvs, allLinxGermlineSvAnnotations, first.svId());
            boolean meetsMaxLength = false;
            if(sv != null)
            {
                meetsMaxLength = Math.abs(sv.start().position() - sv.end().position()) <= MAX_LENGTH_FOR_IMPLIED_DELS;
            }

            boolean hasNoExistingGermlineDel = !hasGermlineDeletionInGene(allGermlineDeletions, first.gene());

            double tumorCopyNumber = Math.max(first.undisruptedCopyNumber(), second.undisruptedCopyNumber());
            GermlineStatus tumorStatus = tumorCopyNumber < 0.5 ? GermlineStatus.HOM_DELETION : GermlineStatus.HET_DELETION;
            boolean wouldBeReportableDeletion = wouldBeReportableGermlineDeletion(first.gene(), tumorStatus, driverGenes);

            if(bothReported && bothDel && sameGene && sameTranscript && meetsMaxLength && hasNoExistingGermlineDel
                    && wouldBeReportableDeletion)
            {
                // assumes deletion is heterozygous in germline
                impliedDeletions.add(new GermlineDeletion(first.gene(), first.chromosome(), first.chromosomeBand(), Math.min(sv.start()
                        .position(), sv.end().position()), Math.max(sv.start().position(), sv.end()
                        .position()), 0, 0, 0, GermlineDetectionMethod.SEGMENT, GermlineStatus.HET_DELETION, tumorStatus, 1D, tumorCopyNumber, Strings.EMPTY, 0, true));
            }
        }
        return impliedDeletions;
    }

    private static boolean wouldBeReportableGermlineDeletion(@NotNull String gene, @NotNull GermlineStatus tumorStatus,
            @NotNull List<DriverGene> driverGenes)
    {
        Optional<DriverGene> optionalDriverGene = driverGenes.stream().filter(d -> d.gene().equals(gene)).findFirst();
        if(optionalDriverGene.isEmpty())
        {
            return false;
        }

        DriverGene driverGene = optionalDriverGene.get();
        if(driverGene.reportGermlineDeletion() == DriverGeneGermlineReporting.NONE)
        {
            return false;
        }
        if(driverGene.reportGermlineDeletion() == DriverGeneGermlineReporting.ANY
                || driverGene.reportGermlineDeletion() == DriverGeneGermlineReporting.VARIANT_NOT_LOST)
        {
            return true;
        }

        if(driverGene.reportGermlineDeletion() == DriverGeneGermlineReporting.WILDTYPE_LOST)
        {
            return tumorStatus == GermlineStatus.HOM_DELETION;
        }

        return false;
    }

    @Nullable
    private static StructuralVariant findBySvId(@NotNull List<StructuralVariant> allPurpleSvs,
            @NotNull List<LinxSvAnnotation> allLinxSvAnnotations, int svIdToFind)
    {
        LinxSvAnnotation match = null;
        int index = 0;
        while(match == null && index < allLinxSvAnnotations.size())
        {
            LinxSvAnnotation svAnnotation = allLinxSvAnnotations.get(index);
            if(svAnnotation.svId() == svIdToFind)
            {
                match = svAnnotation;
            }
            index++;
        }

        if(match != null)
        {
            for(StructuralVariant structuralVariant : allPurpleSvs)
            {
                if(structuralVariant.id().equals(match.vcfId()))
                {
                    return structuralVariant;
                }
            }
        }

        LOGGER.warn("Could not find germline structural variant with svId: {}", svIdToFind);
        return null;
    }

    private static boolean hasGermlineDeletionInGene(@NotNull List<GermlineDeletion> germlineDeletions, @NotNull String geneToFind)
    {
        for(GermlineDeletion deletion : germlineDeletions)
        {
            if(deletion.GeneName.equals(geneToFind))
            {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<PurpleGainLoss> selectReportableGainLosses(@NotNull Map<PurpleGainLoss, Boolean> fullLossToReportability)
    {
        List<PurpleGainLoss> reportable = Lists.newArrayList();
        for(Map.Entry<PurpleGainLoss, Boolean> entry : fullLossToReportability.entrySet())
        {
            PurpleGainLoss gainLoss = entry.getKey();
            boolean reported = entry.getValue();
            if(reported)
            {
                reportable.add(gainLoss);
            }
        }
        return reportable;
    }

    @NotNull
    private static List<PurpleLossOfHeterozygosity> selectReportableLossOfHeterozygosities(
            @NotNull Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability)
    {
        List<PurpleLossOfHeterozygosity> reportable = Lists.newArrayList();
        for(Map.Entry<PurpleLossOfHeterozygosity, Boolean> entry : lossOfHeterozygosityToReportability.entrySet())
        {
            PurpleLossOfHeterozygosity lossOfHeterozygosity = entry.getKey();
            boolean reported = entry.getValue();
            if(reported)
            {
                reportable.add(lossOfHeterozygosity);
            }
        }
        return reportable;
    }

    @NotNull
    private static List<PurpleGainLoss> extractAllGainsLosses(@NotNull Set<PurpleQCStatus> qcStatus, @NotNull Gender gender, double ploidy,
            boolean isTargetRegions, @NotNull List<GeneCopyNumber> allGeneCopyNumbers)
    {
        List<DriverGene> allGenes = Lists.newArrayList();
        for(GeneCopyNumber geneCopyNumber : allGeneCopyNumbers)
        {
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
        DeletionDrivers delDrivers = new DeletionDrivers(qcStatus, allGenesPanel);

        List<DriverCatalog> allGainLosses = Lists.newArrayList();

        allGainLosses.addAll(AmplificationDrivers.findAmplifications(qcStatus, gender, allGenesPanel, ploidy, allGeneCopyNumbers, isTargetRegions));

        allGainLosses.addAll(delDrivers.deletions(allGeneCopyNumbers, isTargetRegions));

        return somaticGainsLossesFromDrivers(allGainLosses);
    }

    @NotNull
    private static List<PurpleGainLoss> somaticGainsLossesFromDrivers(@NotNull List<DriverCatalog> drivers)
    {
        List<PurpleGainLoss> gainsLosses = Lists.newArrayList();

        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(drivers);
        for(DriverCatalogKey key : geneDriverMap.keySet())
        {
            DriverCatalog geneDriver = geneDriverMap.get(key);

            if(geneDriver.driver() == DriverType.AMP || geneDriver.driver() == DriverType.PARTIAL_AMP
                    || geneDriver.driver() == DriverType.DEL)
            {
                gainsLosses.add(toGainLoss(geneDriver));
            }
        }
        return gainsLosses;
    }

    @NotNull
    private static PurpleGainLoss toGainLoss(@NotNull DriverCatalog driver)
    {
        return ImmutablePurpleGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(CopyNumberInterpretationUtil.fromCNADriver(driver))
                .minCopies(Math.max(0, driver.minCopyNumber()))
                .maxCopies(Math.max(0, driver.maxCopyNumber()))
                .build();
    }

    @NotNull
    private static PurpleFit createFit(@NotNull PurpleData purple)
    {
        return ImmutablePurpleFit.builder()
                .qc(PurpleConversion.convert(purple.purityContext().qc()))
                .fittedPurityMethod(PurpleFittedPurityMethod.valueOf(purple.purityContext().method().name()))
                .purity(purple.purityContext().bestFit().purity())
                .minPurity(purple.purityContext().score().minPurity())
                .maxPurity(purple.purityContext().score().maxPurity())
                .ploidy(purple.purityContext().bestFit().ploidy())
                .minPloidy(purple.purityContext().score().minPloidy())
                .maxPloidy(purple.purityContext().score().maxPloidy())
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createCharacteristics(@NotNull PurpleData purple)
    {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(purple.purityContext().wholeGenomeDuplication())
                .microsatelliteIndelsPerMb(purple.purityContext().microsatelliteIndelsPerMb())
                .microsatelliteStatus(PurpleMicrosatelliteStatus.valueOf(purple.purityContext().microsatelliteStatus().name()))
                .tumorMutationalBurdenPerMb(purple.purityContext().tumorMutationalBurdenPerMb())
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext()
                        .tumorMutationalBurdenStatus()
                        .name()))
                .tumorMutationalLoad(purple.purityContext().tumorMutationalLoad())
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext().tumorMutationalLoadStatus().name()))
                .svTumorMutationalBurden(purple.purityContext().svTumorMutationalBurden())
                .build();
    }
}
