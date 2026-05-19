package com.hartwig.hmftools.finding.util;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.DisruptionBuilder;
import com.hartwig.hmftools.finding.datamodel.Doubles;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.FusionBuilder;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class CopyNumberRoundingTransformer
{

    private static final double COPY_NUMBER_ROUNDING_THRESHOLD = 10;
    private static final double COPY_NUMBER_FILTER = 0.1;
    private static final Logger LOGGER = LogManager.getLogger(CopyNumberRoundingTransformer.class);

    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .germlineSmallVariants(roundCopyNumberSmallVariant(record.germlineSmallVariants()))
                .somaticSmallVariants(roundCopyNumberSmallVariant(record.somaticSmallVariants()))
                .somaticDisruptions(roundCopyNumberDisruption(record.somaticDisruptions()))
                .fusions(roundCopyNumberFusion(record.fusions()))
                .hlaAlleles(roundHlACopyNumbers(record.hlaAlleles()))
                .build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> roundCopyNumberSmallVariant(@NotNull DriverFindingList<SmallVariant> variants)
    {
        return TransformUtil.transformDriverFindingList(variants, CopyNumberRoundingTransformer::roundCopyNumbers);
    }

    @NotNull
    private static DriverFindingList<SmallVariant> roundCopyNumberSomaticVariant(@NotNull DriverFindingList<SmallVariant> somaticVariants)
    {
        return DriverFindingListBuilder.builder(somaticVariants).findings(transformSmallVariant(somaticVariants.findings())).build();
    }

    @NotNull
    static SmallVariant roundCopyNumbers(@NotNull SmallVariant smallVariant)
    {
        return SmallVariantBuilder.builder(smallVariant)
                .adjustedCopyNumber(roundCopyNumber(smallVariant.adjustedCopyNumber()))
                .build();
    }

    @NotNull
    private static DriverFindingList<Disruption> roundCopyNumberDisruption(@NotNull DriverFindingList<Disruption> disruptions)
    {
        return TransformUtil.transformDriverFindingList(disruptions, CopyNumberRoundingTransformer::roundCopyNumbers);
    }

    @NotNull
    static Disruption roundCopyNumbers(@NotNull Disruption disruption)
    {
        return filter(DisruptionBuilder.builder(disruption)
                .disruptedCopyNumber(roundCopyNumber(disruption.disruptedCopyNumber()))
                .undisruptedCopyNumber(roundCopyNumber(disruption.undisruptedCopyNumber()))
                .build());
    }

    @NotNull
    static Disruption filter(@NotNull Disruption disruption)
    {
        if(disruption.disruptedCopyNumber() >= COPY_NUMBER_FILTER)
        {
            return disruption;
        }
        else
        {
            LOGGER.warn("Gene disruption {} with type {} has less than < {} copies, investigate further",
                    disruption.gene(),
                    disruption.type(), COPY_NUMBER_FILTER);
            return DisruptionBuilder.builder(disruption)
                    .driver(TransformUtil.toCandidate(disruption.driver()).build())
                    .build();
        }
    }

    @NotNull
    private static DriverFindingList<Fusion> roundCopyNumberFusion(@NotNull DriverFindingList<Fusion> fusions)
    {
        return TransformUtil.transformDriverFindingList(fusions, CopyNumberRoundingTransformer::roundCopyNumbers);
    }

    @NotNull
    static Fusion roundCopyNumbers(@NotNull Fusion fusion)
    {
        return filter(FusionBuilder.builder(fusion)
                .junctionCopyNumber(roundCopyNumber(fusion.junctionCopyNumber()))
                .build());
    }

    @NotNull
    static Fusion filter(@NotNull Fusion fusion)
    {

        if(fusion.junctionCopyNumber() >= COPY_NUMBER_FILTER)
        {
            return fusion;
        }
        else
        {
            String fusionName = fusion.geneUp() + " - " + fusion.geneDown();
            LOGGER.warn("Fusion {} has less than < {} copies, investigate further", fusionName, COPY_NUMBER_FILTER);
            return FusionBuilder.builder(fusion)
                    .driver(TransformUtil.toCandidate(fusion.driver()).build())
                    .build();
        }
    }

    @NotNull
    private static FindingList<HlaAllele> roundHlACopyNumbers(@NotNull FindingList<HlaAllele> hlaAlleles)
    {
        return TransformUtil.transformFindingList(hlaAlleles, CopyNumberRoundingTransformer::roundCopyNumbers);
    }

    @NotNull
    static HlaAllele roundCopyNumbers(HlaAllele hlaAllele) {
        return HlaAlleleBuilder.builder(hlaAllele).tumorCopyNumber(roundCopyNumberNullable(hlaAllele.tumorCopyNumber())).build();
    }

    @Nullable
    static Double roundCopyNumberNullable(@Nullable Double value)
    {
        if(value != null)
        {
            return roundCopyNumber(value);
        }
        return null;
    }

    @NotNull
    static Double roundCopyNumber(@NotNull Double value)
    {
        double doubleValue = Math.max(value, 0);
        return Doubles.round(doubleValue, doubleValue > COPY_NUMBER_ROUNDING_THRESHOLD ? 0 : 1);
    }
}