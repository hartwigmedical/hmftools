package com.hartwig.hmftools.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.finding.datamodel.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.DriverSource;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.finding.datamodel.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.GainDeletionBuilder;
import com.hartwig.hmftools.finding.datamodel.ReportedStatus;

import org.jetbrains.annotations.Nullable;

final class GainDeletionFactory
{

    public static DriverFindingList<GainDeletion> somaticGainDeletionFindings(
            OrangeRefGenomeVersion orangeRefGenomeVersion,
            FindingsStatus findingsStatus,
            PurpleRecord purple)
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
            OrangeRefGenomeVersion orangeRefGenomeVersion,
            PurpleRecord purple)
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
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities =
                Objects.requireNonNull(purple.reportableGermlineLossOfHeterozygosities());
        List<PurpleDriver> germlineDrivers = Objects.requireNonNull(purple.germlineDrivers());

        for(PurpleGainDeletion fullDels : reportableGermlineFullDels)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = Objects.requireNonNull(
                    findDriver(germlineDrivers, fullDels.gene(), fullDels.transcript(), PurpleDriverType.GERMLINE_DELETION));

            final PurpleGeneCopyNumber geneCopyNumber =
                    findPurpleGeneCopyNumber(somaticGeneCopyNumbers, fullDels.gene(), fullDels.transcript());

            gainDeletions.add(toGainDel(fullDels, driver, DriverSource.GERMLINE,
                    GainDeletion.Type.HOM_DEL, GainDeletion.Type.HOM_DEL,
                    geneCopyNumber, cnPerChromosome));
        }

        for(PurpleLossOfHeterozygosity loh : reportableGermlineLossOfHeterozygosities)
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = Objects.requireNonNull(
                    findDriver(germlineDrivers, loh.gene(), loh.transcript(), PurpleDriverType.GERMLINE_DELETION));

            final PurpleGeneCopyNumber geneCopyNumber = findPurpleGeneCopyNumber(somaticGeneCopyNumbers, loh.gene(), loh.transcript());

            gainDeletions.add(germlineLohToGainDel(loh, driver, geneCopyNumber, cnPerChromosome));
        }

        gainDeletions.sort(GainDeletion.COMPARATOR);
        return DriverFindingListBuilder.<GainDeletion>builder()
                .status(FindingsStatus.OK)
                .findings(gainDeletions)
                .build();
    }

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

            final GainDeletion.Type somaticGainDelType = switch(gainDeletion.interpretation())
            {
                case FULL_GAIN, PARTIAL_GAIN -> GainDeletion.Type.GAIN;
                case FULL_DEL, PARTIAL_DEL -> GainDeletion.Type.HOM_DEL;
            };

            final PurpleGeneCopyNumber geneCopyNumber =
                    findPurpleGeneCopyNumber(somaticGeneCopyNumbers, gainDeletion.gene(), gainDeletion.transcript());

            somaticGainsDels.add(toGainDel(gainDeletion, driver, DriverSource.SOMATIC,
                    somaticGainDelType,
                    null,
                    geneCopyNumber,
                    cnPerChromosome));
        }
        return somaticGainsDels;
    }

    public static List<GainDeletion> somaticLoh(List<PurpleGeneCopyNumber> lohGeneCopyNumbers, ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return lohGeneCopyNumbers.stream().map(o -> somaticLohToGainDel(
                o, cnPerChromosome)).toList();
    }

    private static GainDeletion somaticLohToGainDel(PurpleGeneCopyNumber geneCopyNumber, ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(DriverSource.SOMATIC,
                                geneCopyNumber.gene(),
                                PurpleDriverType.DEL,
                                geneCopyNumber.isCanonical(),
                                geneCopyNumber.transcript()))
                        .driverSource(DriverSource.SOMATIC)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.LOW)
                        .driverLikelihood(0)
                        .build()
                )
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.gene())
                .transcript(geneCopyNumber.transcript())
                .isCanonical(geneCopyNumber.isCanonical())
                .somaticType(GainDeletion.Type.HET_DEL)
                .geneExtent(GainDeletion.GeneExtent.FULL_GENE) // not strictly correct
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
            DriverSource sourceSample,
            GainDeletion.Type somaticType,
            @Nullable GainDeletion.Type germlineType,
            PurpleGeneCopyNumber geneCopyNumber,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(sourceSample,
                                purpleGainDeletion.gene(),
                                driver.type(),
                                driver.isCanonical(),
                                purpleGainDeletion.transcript()))
                        .driverSource(sourceSample)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                        .driverLikelihood(driver.driverLikelihood())
                        .build()
                )
                .chromosome(purpleGainDeletion.chromosome())
                .chromosomeBand(purpleGainDeletion.chromosomeBand())
                .gene(purpleGainDeletion.gene())
                .transcript(purpleGainDeletion.transcript())
                .isCanonical(driver.isCanonical())
                .somaticType(somaticType)
                .germlineType(germlineType)
                .geneExtent(toGeneExtent(purpleGainDeletion.interpretation()))
                .tumorMinCopies(purpleGainDeletion.minCopies())
                .tumorMaxCopies(purpleGainDeletion.maxCopies())
                .tumorMinMinorAlleleCopies(geneCopyNumber.minMinorAlleleCopyNumber())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(geneCopyNumber.chromosome(), geneCopyNumber.chromosomeBand()))
                .build();
    }

    private static GainDeletion germlineLohToGainDel(PurpleLossOfHeterozygosity loh, final PurpleDriver driver,
            PurpleGeneCopyNumber geneCopyNumber,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {

        GainDeletion.GeneExtent geneExtent = switch(loh.geneProportion())
        {
            case FULL_GENE -> GainDeletion.GeneExtent.FULL_GENE;
            case PARTIAL_GENE -> GainDeletion.GeneExtent.PARTIAL_GENE;
        };

        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(DriverSource.GERMLINE,
                                loh.gene(),
                                PurpleDriverType.GERMLINE_DELETION,
                                driver.isCanonical(),
                                loh.transcript()))
                        .driverSource(DriverSource.GERMLINE)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                        .driverLikelihood(driver.driverLikelihood())
                        .build()
                )
                .chromosome(loh.chromosome())
                .chromosomeBand(loh.chromosomeBand())
                .gene(loh.gene())
                .transcript(loh.transcript())
                .isCanonical(driver.isCanonical())
                .somaticType(GainDeletion.Type.HET_DEL)
                .germlineType(GainDeletion.Type.HET_DEL)
                .geneExtent(geneExtent)
                .tumorMinCopies(loh.minCopies())
                .tumorMaxCopies(loh.maxCopies())
                .tumorMinMinorAlleleCopies(geneCopyNumber.minMinorAlleleCopyNumber())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(loh.chromosome(), loh.chromosomeBand()))
                .build();
    }

    private static GainDeletion.GeneExtent toGeneExtent(CopyNumberInterpretation copyNumberInterpretation) {
        return switch (copyNumberInterpretation) {
            case FULL_GAIN, FULL_DEL -> GainDeletion.GeneExtent.FULL_GENE;
            case PARTIAL_GAIN, PARTIAL_DEL -> GainDeletion.GeneExtent.PARTIAL_GENE;
        };
    }
}