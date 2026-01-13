package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

class FindingUtil {

    @NotNull
    static <T extends Finding> FindingList<T> notAvailableFindingList() {
        return ImmutableFindingList.<T>builder()
                .status(FindingsStatus.NOT_AVAILABLE)
                .build();
    }

    @NotNull
    static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList() {
        return ImmutableDriverFindingList.<T>builder()
                .status(FindingsStatus.NOT_AVAILABLE)
                .build();
    }

    @NotNull
    static <T> FindingItem<T> notAvailableFindingItem() {
        return ImmutableFindingItem.<T>builder()
                .status(FindingsStatus.NOT_AVAILABLE)
                .build();
    }
}
