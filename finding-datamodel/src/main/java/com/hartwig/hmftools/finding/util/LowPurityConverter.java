package com.hartwig.hmftools.finding.util;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.Driver;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.FindingsStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.ResultStatus;

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
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> convert(@NotNull DriverFindingList<T> driverFindingList, boolean isLowPurity)
    {
        if(shouldConvert(driverFindingList.status(), isLowPurity))
        {
            return DriverFindingListBuilder.<T>builder()
                    .status(convert(driverFindingList.status()))
                    .findings(List.of())
                    .build();
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

    private static FindingsStatus convert(FindingsStatus findingsStatus)
    {
        return FindingsStatusBuilder.builder()
                .status(ResultStatus.NOT_RELIABLE)
                .errors(addLowPurity(findingsStatus.errors()))
                .warnings(removeLowPurity(findingsStatus.warnings()))
                .build();
    }

    private static boolean shouldConvert(FindingsStatus findingsStatus, boolean isLowPurity)
    {
        return findingsStatus.status() == ResultStatus.OK && isLowPurity;
    }

    private static SortedSet<ResultIssue> addLowPurity(SortedSet<ResultIssue> sortedSet)
    {
        SortedSet<ResultIssue> result = new TreeSet<>(sortedSet);
        result.add(ResultIssue.LOW_PURITY);
        return result;
    }

    private static SortedSet<ResultIssue> removeLowPurity(SortedSet<ResultIssue> sortedSet)
    {
        SortedSet<ResultIssue> result = new TreeSet<>(sortedSet);
        result.remove(ResultIssue.LOW_PURITY);
        return result;
    }
}
