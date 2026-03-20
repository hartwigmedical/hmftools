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
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.finding.ResultStatus;

class FindingUtil
{
    static <T extends Driver> DriverFindingList<T> refRequired()
    {
        return notAvailableDriverFindingList(Set.of(ResultIssue.REF_REQUIRED));
    }

    static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList(Set<ResultIssue> errors)
    {
        return DriverFindingListBuilder.<T>builder()
                .status(findingStatus(ResultStatus.NOT_AVAILABLE, errors))
                .findings(List.of())
                .build();
    }

    static <T> FindingItem<T> notAvailableFindingItem(Set<ResultIssue> errors)
    {
        return FindingItemBuilder.<T>builder()
                .status(findingStatus(ResultStatus.NOT_AVAILABLE, errors))
                .build();
    }

    static FindingStatus okStatus()
    {
        return findingStatus(ResultStatus.OK, Set.of());
    }

    static FindingStatus notAvailableStatus(Set<ResultIssue> errors)
    {
        return FindingStatusBuilder.builder()
                .status(ResultStatus.NOT_AVAILABLE)
                .errors(new TreeSet<>(errors))
                .warnings(new TreeSet<>())
                .build();
    }
}
