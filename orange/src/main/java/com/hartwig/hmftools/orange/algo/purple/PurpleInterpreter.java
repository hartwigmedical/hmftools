package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.driver.DriverCatalogFactory.createCopyNumberDriver;
import static com.hartwig.hmftools.common.driver.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.driver.DriverCategory.TSG;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogKey;
import com.hartwig.hmftools.common.driver.DriverCatalogMap;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.ReportableStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.purple.ChromosomalRearrangements;
import com.hartwig.hmftools.datamodel.purple.ImmutableChromosomalRearrangements;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
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
import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    private static final int MAX_LENGTH_FOR_IMPLIED_DELS = 1500;
    private final PurpleVariantFactory purpleVariantFactory;
    private final GermlineGainDeletionFactory germlineGainDelFactory;
    private final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory;
    private final List<DriverGene> driverGenes;
    private final LinxRecord linx;
    private final ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer;

    @Nullable
    private final ChordData chord;
    boolean convertGermlineToSomatic;

    public PurpleInterpreter(final PurpleVariantFactory purpleVariantFactory,
            final GermlineGainDeletionFactory germlineGainDeletionFactory,
            final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory,
            final List<DriverGene> driverGenes, final LinxRecord linx,
            ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer,
            @Nullable final ChordData chord,
            boolean convertGermlineToSomatic)
    {
        this.purpleVariantFactory = purpleVariantFactory;
        this.germlineGainDelFactory = germlineGainDeletionFactory;
        this.germlineLossOfHeterozygosityFactory = germlineLossOfHeterozygosityFactory;
        this.driverGenes = driverGenes;
        this.linx = linx;
        this.chromosomalRearrangementsDeterminer = chromosomalRearrangementsDeterminer;
        this.chord = chord;
        this.convertGermlineToSomatic = convertGermlineToSomatic;
    }

    public PurpleRecord interpret(final PurpleData purple)
    {
        LOGGER.info("Analysing purple data");

        List<PurpleVariant> allSomaticVariants = purpleVariantFactory.fromPurpleVariantContext(purple.allSomaticVariants());

        List<PurpleVariant> reportableSomaticVariants = purpleVariantFactory.fromPurpleVariantContext(purple.reportableSomaticVariants());

        List<PurpleVariant> additionalSuspectSomaticVariants = SomaticVariantSelector.selectInterestingUnreportedVariants(
                allSomaticVariants, reportableSomaticVariants, driverGenes);

        LOGGER.info(" Found an additional {} somatic variants that are potentially interesting",
                additionalSuspectSomaticVariants.size());

        List<PurpleVariant> allGermlineVariants = purpleVariantFactory.fromPurpleVariantContext(purple.allGermlineVariants());
        List<PurpleVariant> reportableGermlineVariants = purpleVariantFactory.fromPurpleVariantContext(purple.reportableGermlineVariants());
        List<PurpleVariant> additionalSuspectGermlineVariants = GermlineVariantSelector.selectInterestingUnreportedVariants(
                allGermlineVariants);

        if(additionalSuspectGermlineVariants != null)
        {
            LOGGER.info(" Found an additional {} germline variants that are potentially interesting",
                    additionalSuspectGermlineVariants.size());
        }

        List<PurpleGainDeletion> allSomaticGainsDels = extractAllGainsDels(purple.allSomaticGeneCopyNumbers());
        List<PurpleGainDeletion> reportableSomaticGainsDels = somaticGainsDelsFromDrivers(purple.somaticDrivers());

        List<PurpleGainDeletion> nearReportableSomaticGains = CopyNumberSelector.selectNearReportableSomaticGains(
                purple.allSomaticGeneCopyNumbers(), purple.purityContext().bestFit().ploidy(), allSomaticGainsDels, driverGenes);
        LOGGER.info(" Found an additional {} near-reportable somatic gains that are potentially interesting",
                nearReportableSomaticGains.size());

        List<PurpleGainDeletion> additionalSuspectSomaticGainsDels = CopyNumberSelector.selectInterestingUnreportedGainsDels(
                allSomaticGainsDels, reportableSomaticGainsDels);
        LOGGER.info(" Found an additional {} somatic gains/deletions that are potentially interesting",
                additionalSuspectSomaticGainsDels.size());

        List<GermlineDeletion> allGermlineDeletions = purple.allGermlineDeletions();
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(purple.allSomaticGeneCopyNumbers(),
                        allGermlineDeletions, purple.purityContext().microsatelliteStatus(), chord != null ? chord.hrStatus() : null);
        LOGGER.info(" Found an additional {} suspect gene copy numbers with LOH", suspectGeneCopyNumbersWithLOH.size());

        List<PurpleGainDeletion> allGermlineFullDels = null;
        List<PurpleGainDeletion> reportableGermlineFullDels = null;
        List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities = null;
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = null;

        if(allGermlineDeletions != null)
        {
            List<GermlineDeletion> impliedDeletions = implyDeletionsFromBreakends(allGermlineDeletions, linx.reportableGermlineBreakends(),
                    purple.allPassingGermlineStructuralVariants(), linx.allGermlineStructuralVariants(), driverGenes);
            LOGGER.info(" Implied {} additional reportable germline deletions from breakends", impliedDeletions.size());

            List<GermlineDeletion> mergedGermlineDeletions = Lists.newArrayList();
            mergedGermlineDeletions.addAll(allGermlineDeletions);
            mergedGermlineDeletions.addAll(impliedDeletions);

            Map<PurpleGainDeletion, Boolean> fullDelToReportability =
                    germlineGainDelFactory.getReportabilityMap(mergedGermlineDeletions, purple.allSomaticGeneCopyNumbers());

            allGermlineFullDels = Lists.newArrayList(fullDelToReportability.keySet());
            reportableGermlineFullDels = selectReportablegainDels(fullDelToReportability);

            LOGGER.info(" Resolved {} germline deletions of which {} are reportable",
                    allGermlineFullDels.size(), reportableGermlineFullDels.size());

            Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability =
                    germlineLossOfHeterozygosityFactory.getReportabilityMap(mergedGermlineDeletions, purple.allSomaticGeneCopyNumbers());

            allGermlineLossOfHeterozygosities = Lists.newArrayList(lossOfHeterozygosityToReportability.keySet());
            reportableGermlineLossOfHeterozygosities = selectReportableLossOfHeterozygosities(lossOfHeterozygosityToReportability);

            LOGGER.info(" Resolved {} germline heterozygous deletions of which {} are reportable",
                    allGermlineLossOfHeterozygosities.size(), reportableGermlineLossOfHeterozygosities.size());
        }

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purple))
                .tumorStats(TumorStatsFactory.compute(purple))
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
                .allSomaticGainsDels(allSomaticGainsDels)
                .reportableSomaticGainsDels(reportableSomaticGainsDels)
                .nearReportableSomaticGains(nearReportableSomaticGains)
                .additionalSuspectSomaticGainsDels(additionalSuspectSomaticGainsDels)
                .allGermlineDeletions(ConversionUtil.mapToIterable(purple.allGermlineDeletions(), PurpleConversion::convert))
                .allGermlineFullDels(allGermlineFullDels)
                .reportableGermlineFullDels(reportableGermlineFullDels)
                .allGermlineLossOfHeterozygosities(allGermlineLossOfHeterozygosities)
                .reportableGermlineLossOfHeterozygosities(reportableGermlineLossOfHeterozygosities)
                .chromosomalRearrangements(createChromosomalRearrangements(purple, chromosomalRearrangementsDeterminer))
                .build();
    }

    @VisibleForTesting
    static List<GermlineDeletion> implyDeletionsFromBreakends(final List<GermlineDeletion> allGermlineDeletions,
            @Nullable List<LinxBreakend> reportableGermlineBreakends, final List<StructuralVariant> allPurpleGermlineSvs,
            @Nullable List<LinxSvAnnotation> allLinxGermlineSvAnnotations, final List<DriverGene> driverGenes)
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
                impliedDeletions.add(new GermlineDeletion(first.gene(), first.chromosome(), first.chromosomeBand(),
                        Math.min(sv.start().position(), sv.end().position()), Math.max(sv.start().position(), sv.end().position()),
                        0, 0, 0, GermlineDetectionMethod.SEGMENT, GermlineStatus.HET_DELETION, tumorStatus, 1D,
                        tumorCopyNumber, Strings.EMPTY, 0, true));
            }
        }
        return impliedDeletions;
    }

    private static boolean wouldBeReportableGermlineDeletion(
            final String gene, final GermlineStatus tumorStatus, final List<DriverGene> driverGenes)
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
    private static StructuralVariant findBySvId(final List<StructuralVariant> allPurpleSvs,
            final List<LinxSvAnnotation> allLinxSvAnnotations, int svIdToFind)
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

    private static boolean hasGermlineDeletionInGene(final List<GermlineDeletion> germlineDeletions, final String geneToFind)
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

    private static List<PurpleGainDeletion> selectReportablegainDels(final Map<PurpleGainDeletion, Boolean> fullDelToReportability)
    {
        List<PurpleGainDeletion> reportable = Lists.newArrayList();
        for(Map.Entry<PurpleGainDeletion, Boolean> entry : fullDelToReportability.entrySet())
        {
            PurpleGainDeletion gainDel = entry.getKey();
            boolean reported = entry.getValue();
            if(reported)
            {
                reportable.add(gainDel);
            }
        }
        return reportable;
    }

    private static List<PurpleLossOfHeterozygosity> selectReportableLossOfHeterozygosities(
            final Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability)
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

    private static List<PurpleGainDeletion> extractAllGainsDels(final List<GeneCopyNumber> allGeneCopyNumbers)
    {
        List<DriverCatalog> allGainDels = Lists.newArrayList();

        for(GeneCopyNumber geneCopyNumber : allGeneCopyNumbers)
        {
            if(geneCopyNumber.reportableStatus() != ReportableStatus.NONE)
            {
                DriverType type = geneCopyNumber.driverType();

                DriverCategory category;
                LikelihoodMethod likelihoodMethod;
                boolean biallelic;
                double likelihood;

                if(type == DriverType.AMP || type == DriverType.PARTIAL_AMP)
                {
                    likelihoodMethod = LikelihoodMethod.AMP;
                    category = ONCO;
                    biallelic = false;
                    likelihood = 1;
                }
                else
                {
                    likelihoodMethod = LikelihoodMethod.DEL;
                    category = TSG;
                    biallelic = type == DriverType.DEL;
                    likelihood = type == DriverType.DEL ? 1 : 0;
                }

                DriverCatalog driverCatalog = createCopyNumberDriver(category, type, likelihoodMethod, biallelic, likelihood, geneCopyNumber);
                allGainDels.add(driverCatalog);
            }
        }

        return somaticGainsDelsFromDrivers(allGainDels);
    }

    private static final Set<DriverType> AMP_DEL_TYPES = Sets.newHashSet(DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL);

    private static List<PurpleGainDeletion> somaticGainsDelsFromDrivers(final List<DriverCatalog> drivers)
    {
        List<PurpleGainDeletion> gainsDels = Lists.newArrayList();

        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(drivers);
        for(DriverCatalogKey key : geneDriverMap.keySet())
        {
            DriverCatalog geneDriver = geneDriverMap.get(key);

            if(AMP_DEL_TYPES.contains(geneDriver.driver()))
            {
                gainsDels.add(toGainDel(geneDriver));
            }
        }
        return gainsDels;
    }

    private static PurpleGainDeletion toGainDel(final DriverCatalog driver)
    {
        return ImmutablePurpleGainDeletion.builder()
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

    private static PurpleFit createFit(final PurpleData purple)
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

    private static PurpleCharacteristics createCharacteristics(final PurpleData purple)
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

    private static ChromosomalRearrangements createChromosomalRearrangements(
            final PurpleData purple, final ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer)
    {
        return ImmutableChromosomalRearrangements.builder()
                .hasTrisomy1q(chromosomalRearrangementsDeterminer.determine1qTrisomy(purple.allSomaticCopyNumbers()))
                .hasCodeletion1p19q(chromosomalRearrangementsDeterminer.determine1p19qCodeletion(purple.allSomaticCopyNumbers()))
                .build();
    }
}
