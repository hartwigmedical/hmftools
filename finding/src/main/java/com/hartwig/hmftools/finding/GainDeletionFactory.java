package com.hartwig.hmftools.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.finding.datamodel.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.DriverSource;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.GermlineAmpDelFields;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.finding.datamodel.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.GainDeletionBuilder;

import org.jspecify.annotations.Nullable;

final class GainDeletionFactory
{
    public static DriverFindingList<GainDeletion> somaticGainDeletionFindings(
            OrangeRefGenomeVersion orangeRefGenomeVersion,
            FindingsStatus findingsStatus,
            PurpleRecord purple)
    {
        List<GainDeletion> gainDeletions = somaticDriverGainDels(purple.somaticGainsDels());

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
            return FindingUtil.notAvailableDriverFindingList();
        }

        List<GainDeletion> gainDeletions = new ArrayList<>();

        for(PurpleGainDeletion gainDels : Objects.requireNonNull(purple.germlineGainsDels()))
        {
            gainDeletions.add(toGainDel(gainDels, DriverSource.GERMLINE));
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

    private static List<GainDeletion> somaticDriverGainDels(List<PurpleGainDeletion> gainDeletions)
    {
        List<GainDeletion> somaticGainsDels = new ArrayList<>();
        for(PurpleGainDeletion gainDeletion : gainDeletions)
        {
            somaticGainsDels.add(toGainDel(gainDeletion, DriverSource.SOMATIC));
        }
        return somaticGainsDels;
    }

    private static GainDeletion toGainDel(PurpleGainDeletion purpleGainDeletion, DriverSource sourceSample)
    {
        PurpleDriver driver = purpleGainDeletion.driver();
        GermlineAmpDelFields germlineAmpDelFields = purpleGainDeletion.germlineAmpDelFields();
        PurpleGermlineStatus somaticStatus = germlineAmpDelFields != null ? germlineAmpDelFields.somaticStatus() : null;
        PurpleGermlineStatus germlineStatus = germlineAmpDelFields != null ? germlineAmpDelFields.germlineStatus() : null;

        GainDeletion.Type somaticType = somaticGainDelType(driver.type(), somaticStatus);
        GainDeletion.Type germlineType = germlineStatus != null ? germlineGainDelType(driver.type(), germlineStatus) : null;

        return GainDeletionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.gainDeletion(sourceSample,
                                purpleGainDeletion.gene(),
                                driver.type(),
                                driver.isCanonical(),
                                purpleGainDeletion.transcript()))
                        .driverSource(sourceSample)
                        .reportedStatus(DriverUtil.reportedStatus(purpleGainDeletion.driver().reportedStatus()))
                        .driverInterpretation(DriverInterpretation.valueOf(driver.driverInterpretation().name()))
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
                .geneExtent(interpretGeneExtent(purpleGainDeletion.geneRange()))
                .exonRange(interpretExonRange(purpleGainDeletion.exonStart(), purpleGainDeletion.exonEnd()))
                .tumorMinCopyNumber(purpleGainDeletion.minCopyNumber())
                .tumorMaxCopyNumber(purpleGainDeletion.maxCopyNumber())
                .tumorMinMinorAlleleCopyNumber(purpleGainDeletion.minMinorAlleleCopies())
                .chromosomeArmCopyNumber(2) // TODO fix later
                .germlineMinCopyNumber(germlineAmpDelFields != null ? germlineAmpDelFields.germlineMinCopyNumber() : null)
                .tpm(purpleGainDeletion.tpm())
                .tpmPercentile(purpleGainDeletion.tpmPercentile())
                .tpmFoldChange(purpleGainDeletion.tpmFoldChange())
                .build();
    }

    private static GainDeletion.GeneExtent interpretGeneExtent(String geneRange)
    {
        return geneRange.equals("FULL") ? GainDeletion.GeneExtent.FULL_GENE : GainDeletion.GeneExtent.PARTIAL_GENE;
    }

    private static GainDeletion.@Nullable ExonRange interpretExonRange(@Nullable Integer exonStart, @Nullable Integer exonEnd)
    {
        if (exonStart == null && exonEnd == null)
        {
            return null;
        }

        return new GainDeletion.ExonRange(exonStart, exonEnd);
    }

    private static GainDeletion.Type germlineGainDelType(PurpleDriverType purpleDriverType, PurpleGermlineStatus purpleGermlineStatus)
    {
        return switch (purpleDriverType) {
            case GERMLINE_DELETION -> switch (purpleGermlineStatus) {
                case HOM_DELETION -> GainDeletion.Type.HOM_DEL;
                case HET_DELETION -> GainDeletion.Type.HET_DEL;
                default -> throw new IllegalStateException("Unexpected germline status: " + purpleGermlineStatus.name());
            };
            case GERMLINE_AMP -> GainDeletion.Type.GAIN;
            default -> throw new IllegalStateException("Unexpected driver type: " + purpleDriverType.name());
        };
    }

    private static GainDeletion.Type somaticGainDelType(PurpleDriverType purpleDriverType, @Nullable PurpleGermlineStatus somaticStatus)
    {
        return switch (purpleDriverType) {
            case AMP, PARTIAL_AMP -> GainDeletion.Type.GAIN;
            case DEL -> GainDeletion.Type.HOM_DEL;
            case HET_DEL -> GainDeletion.Type.HET_DEL;
            case LOH -> GainDeletion.Type.CN_NEUTRAL_LOH;
            case GERMLINE_DELETION, GERMLINE_AMP -> switch (Objects.requireNonNull(somaticStatus)) {
                case HOM_DELETION -> GainDeletion.Type.HOM_DEL;
                case HET_DELETION -> GainDeletion.Type.HET_DEL;
                case AMPLIFICATION -> GainDeletion.Type.GAIN;
                case NOISE, DIPLOID, UNKNOWN -> GainDeletion.Type.NONE;
            };
            default -> throw new IllegalStateException("Unexpected driver type: " + purpleDriverType.name());
        };
    }
}