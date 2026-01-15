package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class FindingUtil {

    @NotNull
    static <T extends Finding> FindingList<T> notAvailableFindingList() {
        return new FindingList<>(FindingsStatus.NOT_AVAILABLE, List.of());
    }

    @NotNull
    static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList() {
        return new DriverFindingList<>(FindingsStatus.NOT_AVAILABLE, List.of());
    }

    @NotNull
    static <T> FindingItem<T> notAvailableFindingItem() {
        return new FindingItem<>(FindingsStatus.NOT_AVAILABLE, null);
    }
}
