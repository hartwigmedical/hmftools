package com.hartwig.hmftools.orange.algo.purple;

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
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
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
    private final PurpleVariantFactory mPurpleVariantFactory;
    private final GermlineGainDeletionFactory mGermlineGainDelFactory;
    private final GermlineLossOfHeterozygosityFactory mGermlineLossOfHeterozygosityFactory;
    private final List<DriverGene> mDriverGenes;
    private final LinxRecord mLinx;

    @Nullable
    private final ChordData mChordData;
    private final boolean mConvertGermlineToSomatic;

    public PurpleInterpreter(
            final PurpleVariantFactory purpleVariantFactory,
            final GermlineGainDeletionFactory germlineGainDeletionFactory,
            final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory,
            final List<DriverGene> driverGenes, final LinxRecord linx,
            @Nullable final ChordData chord,
            boolean convertGermlineToSomatic)
    {
        mPurpleVariantFactory = purpleVariantFactory;
        mGermlineGainDelFactory = germlineGainDeletionFactory;
        mGermlineLossOfHeterozygosityFactory = germlineLossOfHeterozygosityFactory;
        mDriverGenes = driverGenes;
        mLinx = linx;
        mChordData = chord;
        mConvertGermlineToSomatic = convertGermlineToSomatic;
    }

    public PurpleRecord interpret(final PurpleData purple)
    {
        LOGGER.info("Analysing purple data");

        List<PurpleVariant> allSomaticVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.allSomaticVariants());

        List<PurpleVariant> driverSomaticVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.driverSomaticVariants());

        List<PurpleVariant> allGermlineVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.allGermlineVariants());
        List<PurpleVariant> driverGermlineVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.driverGermlineVariants());

        List<PurpleGainDeletion> driverSomaticGainsDels = somaticGainsDelsFromDrivers(purple.somaticDrivers());

        List<GermlineDeletion> allGermlineDeletions = purple.allGermlineDeletions();

        List<PurpleGainDeletion> allGermlineFullDels = null;
        List<PurpleGainDeletion> driverGermlineDeletions = null;
        List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities = null;
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = null;

        if(allGermlineDeletions != null)
        {
            List<GermlineDeletion> impliedDeletions = implyDeletionsFromBreakends(allGermlineDeletions, mLinx.driverGermlineBreakends(),
                    purple.allPassingGermlineStructuralVariants(), mLinx.allGermlineStructuralVariants(), mDriverGenes);
            LOGGER.info(" Implied {} additional reportable germline deletions from breakends", impliedDeletions.size());

            List<GermlineDeletion> mergedGermlineDeletions = Lists.newArrayList();
            mergedGermlineDeletions.addAll(allGermlineDeletions);
            mergedGermlineDeletions.addAll(impliedDeletions);

            Map<PurpleGainDeletion, Boolean> fullDelToReportability =
                    mGermlineGainDelFactory.getReportabilityMap(mergedGermlineDeletions, purple.somaticGeneCopyNumbers());

            allGermlineFullDels = Lists.newArrayList(fullDelToReportability.keySet());
            driverGermlineDeletions = selectReportablegainDels(fullDelToReportability);

            LOGGER.info(" Resolved {} germline deletions of which {} are reportable",
                    allGermlineFullDels.size(), driverGermlineDeletions.size());

            Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability =
                    mGermlineLossOfHeterozygosityFactory.getReportabilityMap(mergedGermlineDeletions, purple.somaticGeneCopyNumbers());

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
                .otherSomaticVariants(allSomaticVariants)
                .driverSomaticVariants(driverSomaticVariants)
                .otherGermlineVariants(allGermlineVariants)
                .driverGermlineVariants(driverGermlineVariants)
                .somaticCopyNumbers(ConversionUtil.mapToIterable(purple.somaticCopyNumbers(), PurpleConversion::convert))
                .somaticGeneCopyNumbers(ConversionUtil.mapToIterable(purple.somaticGeneCopyNumbers(), PurpleConversion::convert))
                .driverSomaticGainsDels(driverSomaticGainsDels)
                .otherGermlineDeletions(ConversionUtil.mapToIterable(purple.allGermlineDeletions(), PurpleConversion::convert))
                .driverGermlineDeletions(driverGermlineDeletions)
                .allGermlineLossOfHeterozygosities(allGermlineLossOfHeterozygosities)
                .driverGermlineLossOfHeterozygosities(reportableGermlineLossOfHeterozygosities)
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

            boolean bothReported = first.reportedStatus() == com.hartwig.hmftools.datamodel.driver.ReportedStatus.REPORTED
                    && second.reportedStatus() == com.hartwig.hmftools.datamodel.driver.ReportedStatus.REPORTED;
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
                        tumorCopyNumber, Strings.EMPTY, 0, ReportedStatus.REPORTED));
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

}
