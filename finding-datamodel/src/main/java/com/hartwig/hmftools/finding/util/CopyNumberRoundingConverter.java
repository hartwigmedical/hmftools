package com.hartwig.hmftools.finding.util;

import java.util.ArrayList;
import java.util.List;

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
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class CopyNumberRoundingConverter
{

    private static final double COPY_NUMBER_ROUNDING_THRESHOLD = 10;

    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .germlineSmallVariants(roundCopyNumberGermlineVariant(record.germlineSmallVariants()))
                .somaticSmallVariants(roundCopyNumberSomaticVariant(record.somaticSmallVariants()))
                .somaticDisruptions(roundCopyNumberDisruption(record.somaticDisruptions()))
                .fusions(roundCopyNumberFusion(record.fusions()))
                .hlaAlleles(roundHlACopyNumbers(record.hlaAlleles()))
                .build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> roundCopyNumberGermlineVariant(@NotNull DriverFindingList<SmallVariant> germlineVariants)
    {
        return DriverFindingListBuilder.builder(germlineVariants).findings(convertSmallVariant(germlineVariants.findings())).build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> roundCopyNumberSomaticVariant(@NotNull DriverFindingList<SmallVariant> somaticVariants)
    {
        return DriverFindingListBuilder.builder(somaticVariants).findings(convertSmallVariant(somaticVariants.findings())).build();
    }

    @NotNull
    public static List<SmallVariant> convertSmallVariant(@NotNull List<SmallVariant> smallVariants)
    {
        List<SmallVariant> smallVariantList = new ArrayList<>();
        for(SmallVariant smallVariant : smallVariants)
        {
            smallVariantList.add(SmallVariantBuilder.builder(smallVariant)
                    .adjustedCopyNumber(roundCopyNumber(smallVariant.adjustedCopyNumber()))
                    .build());
        }
        return smallVariantList;
    }

    @NotNull
    private static DriverFindingList<Disruption> roundCopyNumberDisruption(@NotNull DriverFindingList<Disruption> disruptions)
    {
        return DriverFindingListBuilder.builder(disruptions).findings(convertDisruption(disruptions.findings())).build();
    }

    @NotNull
    public static List<Disruption> convertDisruption(@NotNull List<Disruption> disruptions)
    {
        List<Disruption> disruptionList = new ArrayList<>();

        for(Disruption disruption : disruptions)
        {
            disruptionList.add(DisruptionBuilder.builder(disruption)
                    .disruptedCopyNumber(roundCopyNumber(disruption.disruptedCopyNumber()))
                    .undisruptedCopyNumber(roundCopyNumber(disruption.undisruptedCopyNumber()))
                    .build());

        }
        return disruptionList;
    }

    @NotNull
    private static DriverFindingList<Fusion> roundCopyNumberFusion(@NotNull DriverFindingList<Fusion> fusions)
    {
        return DriverFindingListBuilder.builder(fusions).findings(convertFusion(fusions.findings())).build();
    }

    @NotNull
    public static List<Fusion> convertFusion(@NotNull List<Fusion> fusions)
    {
        List<Fusion> fusionsList = new ArrayList<>();

        for(Fusion fusion : fusions)
        {
            fusionsList.add(FusionBuilder.builder(fusion).junctionCopyNumber(roundCopyNumber(fusion.junctionCopyNumber()))
                    .build());
        }
        return fusionsList;
    }

    @NotNull
    private static FindingList<HlaAllele> roundHlACopyNumbers(@NotNull FindingList<HlaAllele> hlaAlleles)
    {
        return FindingListBuilder.builder(hlaAlleles).findings(convertHla(hlaAlleles.findings())).build();
    }

    @NotNull
    public static List<HlaAllele> convertHla(@NotNull List<HlaAllele> hlaAlleles)
    {
        List<HlaAllele> hlaList = new ArrayList<>();

        for(HlaAllele hlaAllele : hlaAlleles)
        {
            hlaList.add(HlaAlleleBuilder.builder(hlaAllele).tumorCopyNumber(roundCopyNumberNullable(hlaAllele.tumorCopyNumber()))
                    .build());
        }
        return hlaList;
    }

    @NotNull
    public static Double roundCopyNumber(@NotNull Double value)
    {
        double doubleValue = Math.max(value, 0);
        return Doubles.round(doubleValue, doubleValue > COPY_NUMBER_ROUNDING_THRESHOLD ? 0 : 1);
    }

    @Nullable
    public static Double roundCopyNumberNullable(@Nullable Double value)
    {
        if(value != null)
        {
            return roundCopyNumber(value);
        }
        return null;
    }
}
