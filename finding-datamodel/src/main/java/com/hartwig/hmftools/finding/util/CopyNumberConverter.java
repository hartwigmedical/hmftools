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
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class CopyNumberConverter {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberConverter.class);

    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record) {
        return FindingRecordBuilder.builder(record)
                .germlineSmallVariants(filterGermlineVariants(record.germlineSmallVariants()))
                .somaticSmallVariants(filterSomaticVariants(record.somaticSmallVariants()))
                .somaticDisruptions(filterDisruption(record.somaticDisruptions()))
                .fusions(filterFusion(record.fusions()))
                .build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> filterGermlineVariants(@NotNull DriverFindingList<SmallVariant> germlineVariants) {
        return DriverFindingListBuilder.builder(germlineVariants).findings(convertSmallVariant(germlineVariants.findings())).build();
    }

    @NotNull
    private static DriverFindingList<SmallVariant> filterSomaticVariants(@NotNull DriverFindingList<SmallVariant> somaticVariants) {
        return DriverFindingListBuilder.builder(somaticVariants).findings(convertSmallVariant(somaticVariants.findings())).build();
    }

    @NotNull
    public static List<SmallVariant> convertSmallVariant(@NotNull List<SmallVariant> smallVariants) {
        List<SmallVariant> smallVariantList = new ArrayList<>();
        for (SmallVariant smallVariant : smallVariants) {
            Double copyNumber = roundCopyNumber(smallVariant.adjustedCopyNumber());

            if (copyNumber != null && copyNumber >= 0.1) {
                smallVariantList.add(SmallVariantBuilder.builder(smallVariant)
                        .driver(DriverFieldsBuilder.builder(smallVariant.driver()).reportedStatus(ReportedStatus.REPORTED).build())
                        .build());
            } else {
                smallVariantList.add(SmallVariantBuilder.builder(smallVariant)
                        .driver(DriverFieldsBuilder.builder(smallVariant.driver()).reportedStatus(ReportedStatus.CANDIDATE).build())
                        .build());
                LOGGER.warn("Gene {} with variant (coding) {} with variant (protein) {} has lesser than <0.1 copies, investigate further",
                        smallVariant.gene(),
                        smallVariant.transcriptImpact().hgvsCodingImpact(),
                        smallVariant.transcriptImpact().hgvsProteinImpact());
            }
        }
        return smallVariantList;
    }

    @NotNull
    private static DriverFindingList<Disruption> filterDisruption(@NotNull DriverFindingList<Disruption> disruptions) {
        return DriverFindingListBuilder.builder(disruptions).findings(convertDisruption(disruptions.findings())).build();
    }

    @NotNull
    public static List<Disruption> convertDisruption(@NotNull List<Disruption> disruptions) {
        List<Disruption> disruptionList = new ArrayList<>();

        for (Disruption disruption : disruptions) {
            Double copyNumber = roundCopyNumber(disruption.disruptedCopyNumber());
            if (copyNumber != null && copyNumber >= 0.1) {
                disruptionList.add(DisruptionBuilder.builder(disruption)
                        .driver(DriverFieldsBuilder.builder(disruption.driver()).reportedStatus(ReportedStatus.REPORTED).build())
                        .build());
            } else {
                disruptionList.add(DisruptionBuilder.builder(disruption)
                        .driver(DriverFieldsBuilder.builder(disruption.driver()).reportedStatus(ReportedStatus.CANDIDATE).build())
                        .build());
                LOGGER.warn("Gene disruption {} with type {} has lesser than <0.1 copies, investigate further",
                        disruption.gene(),
                        disruption.type());
            }
        }
        return disruptionList;
    }

    @NotNull
    private static DriverFindingList<Fusion> filterFusion(@NotNull DriverFindingList<Fusion> fusions) {
        return DriverFindingListBuilder.builder(fusions).findings(convertFusion(fusions.findings())).build();
    }

    @NotNull
    public static List<Fusion> convertFusion(@NotNull List<Fusion> fusions) {
        List<Fusion> fusionsList = new ArrayList<>();

        for (Fusion fusion : fusions) {
            Double copyNumber = roundCopyNumber(fusion.junctionCopyNumber());
            String fusionName = fusion.geneUp() + " - " + fusion.geneDown();
            if (copyNumber != null && copyNumber >= 0.1) {
                fusionsList.add(FusionBuilder.builder(fusion)
                        .driver(DriverFieldsBuilder.builder(fusion.driver()).reportedStatus(ReportedStatus.REPORTED).build())
                        .build());
            } else {
                fusionsList.add(FusionBuilder.builder(fusion)
                        .driver(DriverFieldsBuilder.builder(fusion.driver()).reportedStatus(ReportedStatus.CANDIDATE).build())
                        .build());
                LOGGER.warn("Fusion {} has lesser than <0.1 copies, investigate further", fusionName);
            }
        }
        return fusionsList;
    }

    //LS: Duplicate of the onco master repo. Is this necessary or can we compare to 0.05?
    @Nullable
    public static Double roundCopyNumber(@Nullable Double value) {
        if (value != null) {
            double doubleValue = Math.max(value, 0);
            return Doubles.round(doubleValue, doubleValue > 10 ? 0 : 1);
        }
        return null;
    }
}
