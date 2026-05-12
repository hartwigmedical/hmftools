package com.hartwig.hmftools.finding.util;

import static java.util.function.Predicate.not;

import java.util.function.Predicate;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;

public class GainDeletionsFilterTransformer
{
    public static FindingRecord transform(FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticGainDeletions(filtered(record.somaticGainDeletions(), not(GainDeletion::isLossOfHeterozygosity)))
                .build();
    }

    private static <T extends Driver> DriverFindingList<T> filtered(DriverFindingList<T> original, Predicate<T> filter)
    {
        return DriverFindingListBuilder.builder(original)
                .findings(original.findings().stream().filter(filter).toList())
                .build();
    }
}
