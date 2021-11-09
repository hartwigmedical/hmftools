package com.hartwig.hmftools.common.doid;

import com.google.common.collect.ListMultimap;

import org.jetbrains.annotations.NotNull;

public final class DoidTestFactory {

    private DoidTestFactory() {
    }

    @NotNull
    public static DoidParents createDoidParents(@NotNull ListMultimap<String, String> relationship) {
        return new DoidParents(relationship);
    }
}
