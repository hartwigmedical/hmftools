package com.hartwig.hmftools.serve.curation;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class DoidLookupTestFactory {

    private DoidLookupTestFactory() {
    }

    @NotNull
    public static DoidLookup dummy() {
        return new DoidLookup(Maps.newHashMap());
    }
}
