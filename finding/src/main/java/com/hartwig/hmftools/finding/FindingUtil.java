package com.hartwig.hmftools.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.finding.Driver;
import com.hartwig.hmftools.datamodel.finding.DriverFindingList;
import com.hartwig.hmftools.datamodel.finding.Finding;
import com.hartwig.hmftools.datamodel.finding.FindingItem;
import com.hartwig.hmftools.datamodel.finding.FindingList;
import com.hartwig.hmftools.datamodel.finding.FindingsStatus;

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
