package com.hartwig.hmftools.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.finding.*;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;

final class GainDeletionFactory
{

    public static DriverFindingList<GainDeletion> somaticGainDeletionFindings(
            OrangeRefGenomeVersion orangeRefGenomeVersion,
            FindingsStatus findingsStatus,
            PurpleRecord purple)
    {
        ChromosomeArmCopyNumberMap cnPerChromosome =
                ChromosomeArmCopyNumberMap.create(purple.somaticCopyNumbers(), orangeRefGenomeVersion);

        List<GainDeletion> gainDeletions = somaticDriverGainDels(purple.somaticGainsDels(), purple.somaticDrivers(),
                purple.somaticGeneCopyNumbers(), cnPerChromosome);

        gainDeletions.sort(GainDeletion.COMPARATOR);

        return DriverFindingListBuilder.<GainDeletion>builder()
                .status(findingsStatus)
                .findings(gainDeletions)
                .build();
    }

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

        List<PurpleGeneCopyNumber> somaticGeneCopyNumbers = purple.somaticGeneCopyNumbers();

        ChromosomeArmCopyNumberMap cnPerChromosome =
                ChromosomeArmCopyNumberMap.create(purple.somaticCopyNumbers(), orangeRefGenomeVersion);

        List<GainDeletion> gainDeletions = new ArrayList<>();
        List<PurpleDriver> germlineDrivers = Objects.requireNonNull(purple.germlineDrivers());

        for(PurpleGainDeletion fullDels : Objects.requireNonNull(purple.germlineGainsDels()))
        {
            // find the purple driver object, it should be there
            PurpleDriver driver = Objects.requireNonNull(
                    findDriver(germlineDrivers, fullDels.gene(), fullDels.transcript(), PurpleDriverType.GERMLINE_DELETION));

            final PurpleGeneCopyNumber geneCopyNumber =
                    findPurpleGeneCopyNumber(somaticGeneCopyNumbers, fullDels.gene(), fullDels.transcript());

            gainDeletions.add(toGainDel(fullDels, driver, GainDeletion.Type.GERMLINE_DEL_HOM_IN_TUMOR, DriverSource.GERMLINE,
                    geneCopyNumber, cnPerChromosome));
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
}