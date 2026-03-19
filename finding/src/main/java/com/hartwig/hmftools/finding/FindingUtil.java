package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.Driver;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.Finding;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingList;
import com.hartwig.hmftools.finding.datamodel.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.FindingsStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.ResultStatus;

class FindingUtil
{
    static <T extends Finding> FindingList<T> emptyFindingList(ResultStatus status)
    {
        return FindingListBuilder.<T>builder()
                .status(findingsStatus(status))
                .findings(List.of())
                .build();
    }

    static <T extends Driver> DriverFindingList<T> emptyDriverFindingList(ResultStatus status)
    {
        return DriverFindingListBuilder.<T>builder()
                .status(findingsStatus(status))
                .findings(List.of())
                .build();
    }

    static <T> FindingItem<T> nullFindingItem(ResultStatus status)
    {
        return FindingItemBuilder.<T>builder()
                .status(findingsStatus(status))
                .build();
    }

    static FindingsStatus findingsStatus(ResultStatus status)
    {
        return FindingsStatusBuilder.builder()
                .status(status)
                .errors(new TreeSet<>())
                .warnings(new TreeSet<>())
                .build();
    }
}
