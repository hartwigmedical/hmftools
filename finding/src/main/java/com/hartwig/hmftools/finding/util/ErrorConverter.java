package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.NO_REPORTABLE_VALUE;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import jakarta.validation.constraints.NotNull;

public class ErrorConverter
{
    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(convert(record.somaticSmallVariants()))
                .germlineSmallVariants(convert(record.germlineSmallVariants()))
                .somaticDisruptions(convert(record.somaticDisruptions()))
                .germlineDisruptions(convert(record.germlineDisruptions()))
                .somaticGainDeletions(convert(record.somaticGainDeletions()))
                .germlineGainDeletions(convert(record.germlineGainDeletions()))
                .fusions(convert(record.fusions()))
                .viruses(convert(record.viruses()))
                .chromosomeArmCopyNumbers(convert(record.chromosomeArmCopyNumbers()))
                // For HLA status remains the same, but tumor fields are cleared.
                .hlaAlleles(convertHla(record.hlaAlleles()))
                .pharmacoGenotypes(convert(record.pharmacoGenotypes()))
                .predictedTumorOrigin(convert(record.predictedTumorOrigin()))
                .microsatelliteStability(convert(record.microsatelliteStability()))
                .tumorMutationalLoad(convert(record.tumorMutationalLoad()))
                .tumorMutationalBurden(convert(record.tumorMutationalBurden()))
                .homologousRecombination(convert(record.homologousRecombination()))
                .build();
    }

    @NotNull
    private static FindingList<HlaAllele> convertHla(@NotNull FindingList<HlaAllele> findingList)
    {
        if(!findingList.status().isOK() || !findingList.findings().isEmpty())
        {
            List<HlaAllele> hlaAlleles = !findingList.status().isOK() ?
                    FindingRecordConverterUtil.convert(findingList.findings(), ErrorConverter::convert, null) : findingList.findings();

            hlaAlleles = hlaAlleles.stream().filter(h -> h.qcStatus().contains(HlaAllele.QcStatus.PASS)).toList();

            FindingStatus findingStatus = findingList.status();
            if(hlaAlleles.isEmpty())
            {
                // Changing status code because this is different from there being no results.
                // The issue is that no results meet the required criteria.
                findingStatus = FindingStatusBuilder.builder(findingStatus)
                        .status(FindingStatus.Status.NOT_AVAILABLE)
                        .errors(new TreeSet<>(Set.of(NO_REPORTABLE_VALUE)))
                        .build();
            }

            return FindingListBuilder.<HlaAllele>builder()
                    .status(findingStatus)
                    .findings(hlaAlleles)
                    .build();
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Finding> FindingList<T> convert(@NotNull FindingList<T> findingList)
    {
        if(!findingList.status().isOK())
        {
            return FindingListBuilder.<T>builder()
                    .status(findingList.status())
                    .findings(List.of())
                    .build();
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> convert(@NotNull DriverFindingList<T> driverFindingList)
    {
        if(!driverFindingList.status().isOK())
        {
            return DriverFindingListBuilder.<T>builder()
                    .status(driverFindingList.status())
                    .findings(List.of())
                    .build();
        }
        else
        {
            return driverFindingList;
        }
    }

    @NotNull
    private static <T> FindingItem<T> convert(@NotNull FindingItem<T> findingItem)
    {
        if(!findingItem.status().isOK())
        {
            return FindingItemBuilder.<T>builder()
                    .status(findingItem.status())
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static HlaAllele convert(HlaAllele hlaAllele)
    {
        return HlaAlleleBuilder.builder(hlaAllele)
                .tumorCopyNumber(null)
                .build();
    }
}
