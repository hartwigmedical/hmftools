package com.hartwig.hmftools.common.doid;

import com.google.common.collect.ListMultimap;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class DoidTestFactory {

    private DoidTestFactory() {
    }

    @NotNull
    public static DoidParents createDoidParents(@NotNull ListMultimap<String, String> relationship) {
        return new DoidParents(relationship);
    }

    @NotNull
    public static DoidNode createDoidNode(@NotNull String doid, @NotNull String doidTerm) {
        return ImmutableDoidNode.builder().doid(doid).url(Strings.EMPTY).doidTerm(doidTerm).build();
    }
}
