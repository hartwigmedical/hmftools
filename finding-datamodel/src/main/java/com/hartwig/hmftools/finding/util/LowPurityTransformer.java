package com.hartwig.hmftools.finding.util;

import java.util.Set;
import java.util.SortedSet;
import java.util.function.Function;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import jakarta.validation.constraints.NotNull;

public class LowPurityTransformer
{
    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        boolean isLowPurity = record.qc().isLowPurity();
        return FindingRecordBuilder.builder(record)
                .somaticDisruptions(transform(record.somaticDisruptions(), isLowPurity))
                .somaticGainDeletions(transform(record.somaticGainDeletions(), isLowPurity))
                .viruses(transform(record.viruses(), isLowPurity))
                .microsatelliteStability(transform(record.microsatelliteStability(), isLowPurity))
                .tumorMutationalLoad(transform(record.tumorMutationalLoad(), isLowPurity))
                .tumorMutationalBurden(transform(record.tumorMutationalBurden(), isLowPurity))
                .homologousRecombination(transform(record.homologousRecombination(), isLowPurity))
                // For HLA status remains the same, but tumor fields are cleared.
                .hlaAlleles(transform(record.hlaAlleles(), isLowPurity, Function.identity(), LowPurityTransformer::transform))
                .build();
    }


    @NotNull
    private static <T extends Finding> FindingList<T> transform(@NotNull FindingList<T> findingList, boolean isLowPurity,
            @NotNull Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @NotNull Function<T, T> findingConverter)
    {
        if(shouldConvert(findingList.status(), isLowPurity))
        {
            return FindingRecordTransformerUtil.transformFindingList(findingList, findingsStatusConverter, findingConverter, null);
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> transform(@NotNull DriverFindingList<T> driverFindingList, boolean isLowPurity)
    {
        if(shouldConvert(driverFindingList.status(), isLowPurity))
        {
            return FindingRecordTransformerUtil.transformDriverFindingList(driverFindingList,
                    LowPurityTransformer::transform,
                    f -> null,
                    null);
        }
        else
        {
            return driverFindingList;
        }
    }

    @NotNull
    private static <T> FindingItem<T> transform(@NotNull FindingItem<T> findingItem, boolean isLowPurity)
    {
        if(shouldConvert(findingItem.status(), isLowPurity))
        {
            return FindingItemBuilder.<T>builder()
                    .status(transform(findingItem.status()))
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static FindingStatus transform(FindingStatus findingStatus)
    {
        return FindingStatusBuilder.builder()
                .status(FindingStatus.Status.NOT_RELIABLE)
                .errors(addLowPurity(findingStatus.errors()))
                .warnings(removeLowPurity(findingStatus.warnings()))
                .build();
    }

    private static boolean shouldConvert(FindingStatus findingStatus, boolean isLowPurity)
    {
        return findingStatus.status() == FindingStatus.Status.OK && isLowPurity;
    }

    private static SortedSet<FindingStatus.Issue> addLowPurity(SortedSet<FindingStatus.Issue> issues)
    {
        return FindingUtil.addIssues(issues, Set.of(FindingStatus.Issue.LOW_PURITY));
    }

    private static SortedSet<FindingStatus.Issue> removeLowPurity(SortedSet<FindingStatus.Issue> issues)
    {
        return FindingUtil.removeIssues(issues, Set.of(FindingStatus.Issue.LOW_PURITY));
    }

    private static HlaAllele transform(HlaAllele hlaAllele)
    {
        return HlaAlleleBuilder.builder(hlaAllele)
                .tumorCopyNumber(null)
                .build();
    }
}
