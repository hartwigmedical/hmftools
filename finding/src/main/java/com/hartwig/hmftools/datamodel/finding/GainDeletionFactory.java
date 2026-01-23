package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;

import org.jetbrains.annotations.NotNull;

final class GainDeletionFactory
{

    public static DriverFindingList<GainDeletion> somaticGainDeletionFindings(
            @NotNull OrangeRefGenomeVersion orangeRefGenomeVersion,
            @NotNull FindingsStatus findingsStatus,
            @NotNull PurpleRecord purple)
    {
        ChromosomeArmCopyNumberMap cnPerChromosome =
                ChromosomeArmCopyNumberMap.create(purple.allSomaticCopyNumbers(), orangeRefGenomeVersion);

        List<GainDeletion> gainDeletions = new ArrayList<>();
        gainDeletions.addAll(somaticDriverGainDels(purple.reportableSomaticGainsDels(), purple.somaticDrivers(), purple.allSomaticGeneCopyNumbers(), cnPerChromosome));

        // we are going to add somatic LOH to purple. For this backported version we will reverse engineer how they might look
        gainDeletions.addAll(somaticLoh(purple.suspectGeneCopyNumbersWithLOH(), cnPerChromosome));

        gainDeletions.sort(GainDeletion.COMPARATOR);

        return DriverFindingListBuilder.<GainDeletion>builder()
                .status(findingsStatus)
                .findings(gainDeletions)
                .build();
    }

    // in orange data, HOM_DELS are stored as germline full dels, HET_DELS are stored in LOH, they do not overlap.
    // findings the reportable ones are in purple drivers. Other types are not reportable, we can ignore them
    public static DriverFindingList<GainDeletion> germlineGainDeletionFindings(
            boolean hasGermlineSample,
            @NotNull OrangeRefGenomeVersion orangeRefGenomeVersion,
            @NotNull PurpleRecord purple)
    {
        if(!hasGermlineSample)
        {
            return DriverFindingListBuilder.<GainDeletion>builder()
                    .status(FindingsStatus.NOT_AVAILABLE)
                    .findings(List.of())
                    .build();
        }

        List<PurpleGeneCopyNumber> somaticGeneCopyNumbers = purple.allSomaticGeneCopyNumbers();

        ChromosomeArmCopyNumberMap cnPerChromosome =
                ChromosomeArmCopyNumberMap.create(purple.allSomaticCopyNumbers(), orangeRefGenomeVersion);

        List<GainDeletion> gainDeletions = new ArrayList<>();
        List<PurpleGainDeletion> reportableGermlineFullDels = Objects.requireNonNull(purple.reportableGermlineFullDels());
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = Objects.requireNonNull(purple.reportableGermlineLossOfHeterozygosities());
        List<PurpleDriver> germlineDrivers = Objects.requireNonNull(purple.germlineDrivers());

        for(PurpleGainDeletion fullDels : reportableGermlineFullDels)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = Objects.requireNonNull(
                    findDriver(germlineDrivers, fullDels.gene(), fullDels.transcript(), PurpleDriverType.GERMLINE_DELETION));

            final PurpleGeneCopyNumber geneCopyNumber =
                    findPurpleGeneCopyNumber(somaticGeneCopyNumbers, fullDels.gene(), fullDels.transcript());

            gainDeletions.add(toGainDel(fullDels, driver, GainDeletion.Type.GERMLINE_DEL_HOM_IN_TUMOR, DriverSource.GERMLINE, geneCopyNumber, cnPerChromosome));
        }

        for(PurpleLossOfHeterozygosity loh : reportableGermlineLossOfHeterozygosities)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = Objects.requireNonNull(
                    findDriver(germlineDrivers, loh.gene(), loh.transcript(), PurpleDriverType.GERMLINE_DELETION));

            final PurpleGeneCopyNumber geneCopyNumber = findPurpleGeneCopyNumber(somaticGeneCopyNumbers, loh.gene(), loh.transcript());

            gainDeletions.add(toGainDel(loh, driver, geneCopyNumber, cnPerChromosome));
        }

        gainDeletions.sort(GainDeletion.COMPARATOR);
        return DriverFindingListBuilder.<GainDeletion>builder()
                .status(FindingsStatus.OK)
                .findings(gainDeletions)
                .build();
    }

    @NotNull
    private static PurpleGeneCopyNumber findPurpleGeneCopyNumber(final List<PurpleGeneCopyNumber> somaticGeneCopyNumbers,
            final String gene, final String transcript)
    {
        return somaticGeneCopyNumbers.stream()
                .filter(o -> o.gene().equals(gene) && o.transcript().equals(transcript))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("No gene copy number found for " + gene + " transcript " + transcript));
    }

    private static List<GainDeletion> somaticDriverGainDels(
            List<PurpleGainDeletion> gainDeletions, final List<PurpleDriver> drivers,
            List<PurpleGeneCopyNumber> somaticGeneCopyNumbers,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        List<GainDeletion> somaticGainsDels = new ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            // we have to reverse the copy number interpretation logic to get back the purple driver type
            final PurpleDriverType purpleDriverType = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN -> PurpleDriverType.AMP;
                case PARTIAL_GAIN -> PurpleDriverType.PARTIAL_AMP;
                case FULL_DEL, PARTIAL_DEL -> PurpleDriverType.DEL;
            };

            PurpleDriver driver = findDriver(drivers, gainDeletion.gene(), gainDeletion.transcript(), purpleDriverType);

            final GainDeletion.Type type = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN, PARTIAL_GAIN -> GainDeletion.Type.SOMATIC_GAIN;
                case FULL_DEL, PARTIAL_DEL -> GainDeletion.Type.SOMATIC_DEL;
            };

            final PurpleGeneCopyNumber geneCopyNumber =
                    findPurpleGeneCopyNumber(somaticGeneCopyNumbers, gainDeletion.gene(), gainDeletion.transcript());

            somaticGainsDels.add(toGainDel(gainDeletion, driver, type, DriverSource.SOMATIC, geneCopyNumber, cnPerChromosome));
        }
        return somaticGainsDels;
    }

    public static List<GainDeletion> somaticLoh(List<PurpleGeneCopyNumber> lohGeneCopyNumbers, ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return lohGeneCopyNumbers.stream().map(o -> toGainDel(o, cnPerChromosome)).toList();
    }

    private static GainDeletion toGainDel(PurpleGeneCopyNumber geneCopyNumber, ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(DriverSource.SOMATIC,
                                geneCopyNumber.gene(),
                                CopyNumberInterpretation.FULL_DEL,
                                geneCopyNumber.isCanonical(),
                                geneCopyNumber.transcript()))
                        .driverSource(DriverSource.SOMATIC)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.LOW)
                        .driverLikelihood(0)
                        .build()
                )
                .type(GainDeletion.Type.SOMATIC_LOH)
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.gene())
                .transcript(geneCopyNumber.transcript())
                .isCanonical(geneCopyNumber.isCanonical())
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .tumorMinCopies(geneCopyNumber.minCopyNumber())
                .tumorMaxCopies(geneCopyNumber.maxCopyNumber())
                .tumorMinMinorAlleleCopies(geneCopyNumber.minMinorAlleleCopyNumber())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(geneCopyNumber.chromosome(), geneCopyNumber.chromosomeBand()))
                .build();
    }

    private static PurpleDriver findDriver(final List<PurpleDriver> drivers, final String gene, final String transcript,
            final PurpleDriverType purpleDriverType)
    {
        return drivers.stream()
                .filter(o -> o.gene().equals(gene) && o.transcript().equals(transcript) && o.type().equals(purpleDriverType))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException(
                        "No driver found for " + gene + " transcript " + transcript + " type " + purpleDriverType));
    }

    private static GainDeletion toGainDel(PurpleGainDeletion purpleGainDeletion,
            final PurpleDriver driver,
            GainDeletion.Type type,
            DriverSource sourceSample,
            PurpleGeneCopyNumber geneCopyNumber,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(sourceSample,
                                purpleGainDeletion.gene(),
                                purpleGainDeletion.interpretation(),
                                driver.isCanonical(),
                                purpleGainDeletion.transcript()))
                        .driverSource(sourceSample)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                        .driverLikelihood(driver.driverLikelihood())
                        .build()
                )
                .type(type)
                .chromosome(purpleGainDeletion.chromosome())
                .chromosomeBand(purpleGainDeletion.chromosomeBand())
                .gene(purpleGainDeletion.gene())
                .transcript(purpleGainDeletion.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(purpleGainDeletion.interpretation())
                .tumorMinCopies(purpleGainDeletion.minCopies())
                .tumorMaxCopies(purpleGainDeletion.maxCopies())
                .tumorMinMinorAlleleCopies(geneCopyNumber.minMinorAlleleCopyNumber())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(geneCopyNumber.chromosome(), geneCopyNumber.chromosomeBand()))
                .build();
    }

    private static GainDeletion toGainDel(PurpleLossOfHeterozygosity loh, final PurpleDriver driver, PurpleGeneCopyNumber geneCopyNumber,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {

        CopyNumberInterpretation copyNumberInterpretation = switch(loh.geneProportion())
        {
            case FULL_GENE -> CopyNumberInterpretation.FULL_GAIN;
            case PARTIAL_GENE -> CopyNumberInterpretation.PARTIAL_DEL;
        };

        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(DriverSource.GERMLINE,
                                loh.gene(),
                                copyNumberInterpretation,
                                driver.isCanonical(),
                                loh.transcript()))
                        .driverSource(DriverSource.GERMLINE)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                        .driverLikelihood(driver.driverLikelihood())
                        .build()
                )
                .type(GainDeletion.Type.GERMLINE_DEL_HET_IN_TUMOR)
                .chromosome(loh.chromosome())
                .chromosomeBand(loh.chromosomeBand())
                .gene(loh.gene())
                .transcript(loh.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(copyNumberInterpretation)
                .tumorMinCopies(loh.minCopies())
                .tumorMaxCopies(loh.maxCopies())
                .tumorMinMinorAlleleCopies(geneCopyNumber.minMinorAlleleCopyNumber())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(loh.chromosome(), loh.chromosomeBand()))
                .build();
    }
}