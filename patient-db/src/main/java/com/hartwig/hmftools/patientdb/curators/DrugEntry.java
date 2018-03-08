package com.hartwig.hmftools.patientdb.curators;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class DrugEntry {

    @NotNull
    private final List<String> names;
    @NotNull
    private final String type;
    @NotNull
    private final String canonicalName;

    public DrugEntry(@NotNull final List<String> names, @NotNull final String type, @NotNull final String canonicalName) {
        this.names = names;
        this.type = type;
        this.canonicalName = canonicalName;
    }
}
