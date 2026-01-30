package com.hartwig.hmftools.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.finding.*;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
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

        List<GainDeletion> gainDeletions = somaticDriverGainDels(purple.somaticGainsDels(), cnPerChromosome);

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

        ChromosomeArmCopyNumberMap cnPerChromosome =
                ChromosomeArmCopyNumberMap.create(purple.somaticCopyNumbers(), orangeRefGenomeVersion);

        List<GainDeletion> gainDeletions = new ArrayList<>();

        for(PurpleGainDeletion gainDels : Objects.requireNonNull(purple.germlineGainsDels()))
        {
            gainDeletions.add(toGainDel(gainDels, DriverSource.GERMLINE, cnPerChromosome));
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
            List<PurpleGainDeletion> gainDeletions, ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        List<GainDeletion> somaticGainsDels = new ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            somaticGainsDels.add(toGainDel(gainDeletion, DriverSource.SOMATIC, cnPerChromosome));
        }
        return somaticGainsDels;
    }

    private static GainDeletion toGainDel(PurpleGainDeletion purpleGainDeletion,
            DriverSource sourceSample,
            ChromosomeArmCopyNumberMap cnPerChromosome)
    {
        PurpleDriver driver = purpleGainDeletion.driver();

        GainDeletion.Type type = switch(purpleGainDeletion.driver().type())
        {
            case AMP, PARTIAL_AMP -> GainDeletion.Type.SOMATIC_GAIN;
            case DEL, HET_DEL -> GainDeletion.Type.SOMATIC_DEL;
            case LOH -> GainDeletion.Type.SOMATIC_LOH;
            case GERMLINE_DELETION -> {
                if(purpleGainDeletion.germlineAmpDelFields().somaticStatus() == PurpleGermlineStatus.HOM_DELETION)
                {
                    yield GainDeletion.Type.GERMLINE_DEL_HOM_IN_TUMOR;
                }
                else if(purpleGainDeletion.germlineAmpDelFields().somaticStatus() == PurpleGermlineStatus.HET_DELETION)
                {
                    yield GainDeletion.Type.GERMLINE_DEL_HET_IN_TUMOR;
                }
                else
                {
                    throw new IllegalStateException("Unexpected germline status: " + purpleGainDeletion.germlineAmpDelFields().somaticStatus());
                }
            }
            case GERMLINE_AMP -> GainDeletion.Type.GERMLINE_GAIN;
            default -> throw new IllegalStateException("Unexpected driver type: " + purpleGainDeletion.driver().type());
        };

        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(sourceSample,
                                purpleGainDeletion.gene(),
                                purpleGainDeletion.interpretation(),
                                driver.isCanonical(),
                                purpleGainDeletion.transcript()))
                        .driverSource(sourceSample)
                        .reportedStatus(DriverUtil.reportedStatus(purpleGainDeletion.driver().reportedStatus()))
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
                .tumorMinMinorAlleleCopies(purpleGainDeletion.minMinorAlleleCopies())
                .chromosomeArmCopies(cnPerChromosome.chromosomeArmCopyNumber(
                        purpleGainDeletion.chromosome(), purpleGainDeletion.chromosomeBand()))
                .build();
    }
}