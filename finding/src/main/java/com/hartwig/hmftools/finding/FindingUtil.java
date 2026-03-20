package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.Driver;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.FindingsStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.ResultStatus;

class FindingUtil
{
    static <T extends Driver> DriverFindingList<T> refRequired()
    {
        return notAvailableDriverFindingList(Set.of(ResultIssue.REF_REQUIRED));
    }

    static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList(Set<ResultIssue> errors)
    {
        return DriverFindingListBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .findings(List.of())
                .build();
    }

    static <T> FindingItem<T> notAvailableFindingItem(Set<ResultIssue> errors)
    {
        return FindingItemBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .build();
    }

    static FindingsStatus notAvailableStatus(Set<ResultIssue> errors)
    {
        return FindingsStatusBuilder.builder()
                .status(ResultStatus.NOT_AVAILABLE)
                .errors(new TreeSet<>(errors))
                .warnings(new TreeSet<>())
                .build();
    }
}
