package com.hartwig.hmftools.finding.util;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.Driver;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;

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
                    .status(FindingsStatus.NOT_RELIABLE)
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
                    .status(FindingsStatus.NOT_RELIABLE)
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static boolean shouldConvert(FindingsStatus findingsStatus, boolean isLowPurity)
    {
        return findingsStatus == FindingsStatus.OK && isLowPurity;
    }
}
