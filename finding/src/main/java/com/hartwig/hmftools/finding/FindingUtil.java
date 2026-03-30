package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

public class FindingUtil
{
    static <T extends Driver> DriverFindingList<T> refRequired()
    {
        return notAvailableDriverFindingList(Set.of(FindingStatus.Issue.REF_REQUIRED));
    }

    public static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList(Set<FindingStatus.Issue> errors)
    {
        return DriverFindingListBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .findings(List.of())
                .build();
    }

    static <T> FindingItem<T> notAvailableFindingItem(Set<FindingStatus.Issue> errors)
    {
        return FindingItemBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .build();
    }

    static FindingStatus notAvailableStatus(Set<FindingStatus.Issue> errors)
    {
        return FindingStatusBuilder.builder()
                .status(FindingStatus.Status.NOT_AVAILABLE)
                .errors(new TreeSet<>(errors))
                .warnings(new TreeSet<>())
                .build();
    }
}
