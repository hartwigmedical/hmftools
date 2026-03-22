package com.hartwig.hmftools.finding.util;

import java.util.SortedSet;
import java.util.TreeSet;
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

public class LowPurityConverter
{
    public static FindingRecord convert(FindingRecord record)
    {
        boolean isLowPurity = record.qc().isLowPurity();
        return FindingRecordBuilder.builder(record)
                .somaticDisruptions(convert(record.somaticDisruptions(), isLowPurity))
                .somaticGainDeletions(convert(record.somaticGainDeletions(), isLowPurity))
                .viruses(convert(record.viruses(), isLowPurity))
                .microsatelliteStability(convert(record.microsatelliteStability(), isLowPurity))
                .tumorMutationalLoad(convert(record.tumorMutationalLoad(), isLowPurity))
                .tumorMutationalBurden(convert(record.tumorMutationalBurden(), isLowPurity))
                .homologousRecombination(convert(record.homologousRecombination(), isLowPurity))
                // For HLA status remains the same, but tumor fields are cleared.
                .hlaAlleles(convert(record.hlaAlleles(), isLowPurity, Function.identity(), LowPurityConverter::convert))
                .build();
    }

    @NotNull
    private static <T extends Finding> FindingList<T> convert(@NotNull FindingList<T> findingList, boolean isLowPurity,
            @NotNull Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @NotNull Function<T, T> findingConverter)
    {
        if(shouldConvert(findingList.status(), isLowPurity))
        {
            return FindingsConverter.convert(findingList, findingsStatusConverter, findingConverter, null);
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> convert(@NotNull DriverFindingList<T> driverFindingList, boolean isLowPurity)
    {
        if(shouldConvert(driverFindingList.status(), isLowPurity))
        {
            return FindingsConverter.convert(driverFindingList,
                    LowPurityConverter::convert,
                    f -> null,
                    null);
        }
        else
        {
            return driverFindingList;
        }
    }

    @NotNull
    private static <T> FindingItem<T> convert(@NotNull FindingItem<T> findingItem, boolean isLowPurity)
    {
        if(shouldConvert(findingItem.status(), isLowPurity))
        {
            return FindingItemBuilder.<T>builder()
                    .status(convert(findingItem.status()))
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static FindingStatus convert(FindingStatus findingStatus)
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

    private static SortedSet<FindingStatus.Issue> addLowPurity(SortedSet<FindingStatus.Issue> sortedSet)
    {
        SortedSet<FindingStatus.Issue> result = new TreeSet<>(sortedSet);
        result.add(FindingStatus.Issue.LOW_PURITY);
        return result;
    }

    private static SortedSet<FindingStatus.Issue> removeLowPurity(SortedSet<FindingStatus.Issue> sortedSet)
    {
        SortedSet<FindingStatus.Issue> result = new TreeSet<>(sortedSet);
        result.remove(FindingStatus.Issue.LOW_PURITY);
        return result;
    }

    private static HlaAllele convert(HlaAllele hlaAllele)
    {
        return HlaAlleleBuilder.builder(hlaAllele)
                .tumorCopyNumber(null)
                .build();
    }
}
