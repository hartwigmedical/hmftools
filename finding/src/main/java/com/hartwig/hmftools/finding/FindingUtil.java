package com.hartwig.hmftools.finding;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.Driver;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.Finding;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingList;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;

class FindingUtil
{
    static <T extends Finding> FindingList<T> emptyFindingList(FindingsStatus status)
    {
        return new FindingList<>(status, List.of());
    }

    static <T extends Driver> DriverFindingList<T> emptyDriverFindingList(FindingsStatus status)
    {
        return new DriverFindingList<>(status, List.of());
    }

    static <T> FindingItem<T> nullFindingItem(FindingsStatus status)
    {
        return new FindingItem<>(status, null);
    }
}
