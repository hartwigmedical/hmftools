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
    static <T extends Finding> FindingList<T> notAvailableFindingList()
    {
        return new FindingList<>(FindingsStatus.NOT_AVAILABLE, List.of());
    }

    static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList()
    {
        return new DriverFindingList<>(FindingsStatus.NOT_AVAILABLE, List.of());
    }

    static <T> FindingItem<T> notAvailableFindingItem()
    {
        return new FindingItem<>(FindingsStatus.NOT_AVAILABLE, null);
    }
}
