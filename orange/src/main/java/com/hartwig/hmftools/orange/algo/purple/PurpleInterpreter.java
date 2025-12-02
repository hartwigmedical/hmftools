package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.util.DriverUtils.convertReportedStatus;
import static com.hartwig.hmftools.orange.conversion.PurpleConversion.convert;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableGainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableMicrosatelliteStability;
import com.hartwig.hmftools.datamodel.finding.ImmutableTumorMutationStatus;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;
import com.hartwig.hmftools.orange.algo.util.FindingKeys;

import org.jetbrains.annotations.Nullable;

public class PurpleInterpreter
{
    private final PurpleVariantFactory mPurpleVariantFactory;
    private final GermlineGainDeletionFactory mGermlineGainDelFactory;
    private final GermlineLossOfHeterozygosityFactory mGermlineLossOfHeterozygosityFactory;

    public PurpleInterpreter(
            final PurpleVariantFactory purpleVariantFactory,
            final GermlineGainDeletionFactory germlineGainDeletionFactory,
            final GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory)
    {
        mPurpleVariantFactory = purpleVariantFactory;
        mGermlineGainDelFactory = germlineGainDeletionFactory;
        mGermlineLossOfHeterozygosityFactory = germlineLossOfHeterozygosityFactory;
    }

    public PurpleRecord interpret(final PurpleData purple)
    {
        LOGGER.info("Analysing purple data");

        List<PurpleVariant> allSomaticVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.allSomaticVariants());

        List<PurpleVariant> driverSomaticVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.driverSomaticVariants());

        List<PurpleVariant> allGermlineVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.allGermlineVariants());
        @Nullable List<PurpleVariant> driverGermlineVariants = mPurpleVariantFactory.fromPurpleVariantContext(purple.driverGermlineVariants());

        List<GainDeletion> driverSomaticGainDels = somaticGainsDelsFromDrivers(purple.somaticDrivers());

        List<GermlineDeletion> allGermlineDeletions = purple.allGermlineDeletions();

        List<GainDeletion> allGermlineFullDels = null;
        List<GainDeletion> driverGermlineDeletions = null;
        List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities = null;
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = null;

        if(allGermlineDeletions != null)
        {
            Map<GainDeletion, Boolean> fullDelToReportability = mGermlineGainDelFactory.getReportabilityMap(
                    allGermlineDeletions, purple.somaticGeneCopyNumbers());

            allGermlineFullDels = Lists.newArrayList(fullDelToReportability.keySet());
            driverGermlineDeletions = selectReportableGainDels(fullDelToReportability);

            LOGGER.info(" Resolved {} germline deletions of which {} are reportable",
                    allGermlineFullDels.size(), driverGermlineDeletions.size());

            Map<PurpleLossOfHeterozygosity, Boolean> lossOfHeterozygosityToReportability =
                    mGermlineLossOfHeterozygosityFactory.getReportabilityMap(allGermlineDeletions, purple.somaticGeneCopyNumbers());

            allGermlineLossOfHeterozygosities = Lists.newArrayList(lossOfHeterozygosityToReportability.keySet());
            reportableGermlineLossOfHeterozygosities = selectReportableLossOfHeterozygosities(lossOfHeterozygosityToReportability);

            LOGGER.info(" Resolved {} germline heterozygous deletions of which {} are reportable",
                    allGermlineLossOfHeterozygosities.size(), reportableGermlineLossOfHeterozygosities.size());
        }

        List<PurpleDriver> somaticDrivers = ConversionUtil.mapToList(purple.somaticDrivers(), PurpleConversion::convert);
        List<SmallVariant> somaticSmallVariants = SmallVariantFactory.create(driverSomaticVariants, somaticDrivers);
        List<PurpleDriver> germlineDrivers = ConversionUtil.mapToNullableList(purple.germlineDrivers(), PurpleConversion::convert);
        List<SmallVariant> germlineSmallVariants = null;

        if(germlineDrivers != null)
        {
            germlineSmallVariants = SmallVariantFactory.create(driverGermlineVariants, germlineDrivers);
        }

        return ImmutablePurpleRecord.builder()
                .fit(createFit(purple))
                .tumorStats(TumorStatsFactory.compute(purple))
                .characteristics(createCharacteristics(purple))
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .driverSomaticVariants(driverSomaticVariants)
                .otherSomaticVariants(allSomaticVariants)
                .driverGermlineVariants(driverGermlineVariants)
                .otherGermlineVariants(allGermlineVariants)
                .driverGermlineVariants(driverGermlineVariants)
                .driverSomaticSmallVariants(somaticSmallVariants)
                .driverGermlineSmallVariants(germlineSmallVariants)
                .somaticCopyNumbers(ConversionUtil.mapToIterable(purple.somaticCopyNumbers(), PurpleConversion::convert))
                .somaticGeneCopyNumbers(ConversionUtil.mapToIterable(purple.somaticGeneCopyNumbers(), PurpleConversion::convert))
                .driverSomaticGainsDels(driverSomaticGainDels)
                .otherGermlineDeletions(ConversionUtil.mapToNullableList(purple.allGermlineDeletions(), PurpleConversion::convert))
                .driverGermlineDeletions(driverGermlineDeletions)
                .allGermlineLossOfHeterozygosities(allGermlineLossOfHeterozygosities)
                .driverGermlineLossOfHeterozygosities(reportableGermlineLossOfHeterozygosities)
                .build();
    }

    private static List<GainDeletion> selectReportableGainDels(final Map<GainDeletion, Boolean> fullDelToReportability)
    {
        List<GainDeletion> reportable = Lists.newArrayList();
        for(Map.Entry<GainDeletion, Boolean> entry : fullDelToReportability.entrySet())
        {
            GainDeletion gainDel = entry.getKey();
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

    private static List<GainDeletion> somaticGainsDelsFromDrivers(final List<DriverCatalog> drivers)
    {
        return gainsDelsFromDrivers(drivers, Set.of(DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL));
    }

    private static List<GainDeletion> germlineGainsDelsFromDrivers(final List<DriverCatalog> drivers)
    {
        return gainsDelsFromDrivers(drivers, Set.of(DriverType.GERMLINE_DELETION));
    }

    private static List<GainDeletion> gainsDelsFromDrivers(final List<DriverCatalog> drivers, final Set<DriverType> driverTypes)
    {
        return drivers.stream()
                .filter(o -> driverTypes.contains(o.driver()))
                .map(PurpleInterpreter::toGainDel)
                .toList();
    }

    private static GainDeletion toGainDel(final DriverCatalog driver)
    {
        CopyNumberInterpretation copyNumberInterpretation = CopyNumberInterpretationUtil.fromCNADriver(driver);

        // TODO: transcript fix
        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.gainDeletion(driver.gene(), copyNumberInterpretation, driver.isCanonical(), driver.transcript()))
                .reportedStatus(convertReportedStatus(driver.reportedStatus()))
                .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                .driver(convert(driver))
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(copyNumberInterpretation)
                .minCopies(Math.max(0, driver.minCopyNumber()))
                .maxCopies(Math.max(0, driver.maxCopyNumber()))
                .build();
    }

    private static PurpleFit createFit(final PurpleData purple)
    {
        return ImmutablePurpleFit.builder()
                .qc(convert(purple.purityContext().qc()))
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
                .microsatelliteStability(ImmutableMicrosatelliteStability.builder()
                        .findingKey(FindingKeys.microsatelliteStability(purple.purityContext().microsatelliteStatus()))
                        .microsatelliteIndelsPerMb(purple.purityContext().microsatelliteIndelsPerMb())
                        .microsatelliteStatus(PurpleMicrosatelliteStatus.valueOf(purple.purityContext().microsatelliteStatus().name()))
                        .build())
                .tumorMutationStatus(ImmutableTumorMutationStatus.builder()
                        .findingKey(FindingKeys.tumorMutationStatus())
                        .tumorMutationalBurdenPerMb(purple.purityContext().tumorMutationalBurdenPerMb())
                        .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext()
                                .tumorMutationalBurdenStatus()
                                .name()))
                        .tumorMutationalLoad(purple.purityContext().tumorMutationalLoad())
                        .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.valueOf(purple.purityContext().tumorMutationalLoadStatus().name()))
                        .svTumorMutationalBurden(purple.purityContext().svTumorMutationalBurden())
                        .build()
                )
                .build();
    }

}
