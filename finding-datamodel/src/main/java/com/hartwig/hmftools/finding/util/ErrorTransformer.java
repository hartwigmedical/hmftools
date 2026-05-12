package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.NO_TUMOR;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;

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

public class ErrorTransformer
{
    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(transform(record.somaticSmallVariants()))
                .germlineSmallVariants(transform(record.germlineSmallVariants()))
                .somaticDisruptions(transform(record.somaticDisruptions()))
                .germlineDisruptions(transform(record.germlineDisruptions()))
                .somaticGainDeletions(transform(record.somaticGainDeletions()))
                .germlineGainDeletions(transform(record.germlineGainDeletions()))
                .fusions(transform(record.fusions()))
                .viruses(transform(record.viruses()))
                .chromosomeArmCopyNumbers(transform(record.chromosomeArmCopyNumbers()))
                // For HLA status remains the same, but tumor fields are cleared.
                .hlaAlleles(transformHla(record.hlaAlleles()))
                .pharmacoGenotypes(transform(record.pharmacoGenotypes()))
                .predictedTumorOrigin(transform(record.predictedTumorOrigin()))
                .microsatelliteStability(transform(record.microsatelliteStability()))
                .tumorMutationalLoad(transform(record.tumorMutationalLoad()))
                .tumorMutationalBurden(transform(record.tumorMutationalBurden()))
                .homologousRecombination(transform(record.homologousRecombination()))
                .build();
    }

    @NotNull
    private static FindingList<HlaAllele> transformHla(@NotNull FindingList<HlaAllele> findingList)
    {
        if(!findingList.status().isOK() || !findingList.findings().isEmpty())
        {
            FindingStatus findingStatus = findingList.status();
            List<HlaAllele> hlaAlleles;
            if(!findingStatus.isOK() && findingStatus.errors().contains(NO_TUMOR) )
            {
                // Clear tumor fields if there is no tumor.
                hlaAlleles = FindingRecordTransformerUtil.transform(findingList.findings(), ErrorTransformer::transform, null);
                // Toggle no tumor from error to warning.
                SortedSet<FindingStatus.Issue> errors = FindingUtil.removeIssues(findingStatus.errors(), Set.of(NO_TUMOR));
                findingStatus = FindingStatusBuilder.builder(findingStatus)
                        .status(errors.isEmpty() ? FindingStatus.Status.OK : findingStatus.status())
                        .errors(errors)
                        .warnings(FindingUtil.addIssues(findingStatus.warnings(), Set.of(NO_TUMOR)))
                        .build();
            }
            else
            {
                hlaAlleles = findingList.findings();
            }

            hlaAlleles = hlaAlleles.stream().filter(h -> h.qcStatus().contains(HlaAllele.QcStatus.PASS)).toList();

            if(hlaAlleles.isEmpty())
            {
                // Changing status code because this is different from there being no results.
                // The issue is that no results meet the required criteria.
                findingStatus = FindingUtil.noReportableValueStatus(findingStatus);
            }

            return FindingListBuilder.builder(findingList)
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
    private static <T extends Finding> FindingList<T> transform(@NotNull FindingList<T> findingList)
    {
        if(!findingList.status().isOK())
        {
            return FindingListBuilder.builder(findingList)
                    .findings(List.of())
                    .build();
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> transform(@NotNull DriverFindingList<T> driverFindingList)
    {
        if(!driverFindingList.status().isOK())
        {
            return DriverFindingListBuilder.builder(driverFindingList)
                    .findings(List.of())
                    .build();
        }
        else
        {
            return driverFindingList;
        }
    }

    @NotNull
    private static <T> FindingItem<T> transform(@NotNull FindingItem<T> findingItem)
    {
        if(!findingItem.status().isOK())
        {
            return FindingItemBuilder.builder(findingItem)
                    .finding(null)
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static HlaAllele transform(HlaAllele hlaAllele)
    {
        return HlaAlleleBuilder.builder(hlaAllele)
                .tumorCopyNumber(null)
                .build();
    }
}
