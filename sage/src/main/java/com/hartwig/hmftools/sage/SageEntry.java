package com.hartwig.hmftools.sage;

import java.util.List;

import com.hartwig.hmftools.sage.context.AltContext;

import org.jetbrains.annotations.NotNull;

public class SageEntry {

    private final AltContext normal;
    private final List<AltContext> tumorAltContexts;

    SageEntry(final AltContext normal, final List<AltContext> tumorAltContexts) {
        assert (!tumorAltContexts.isEmpty());

        this.normal = normal;
        this.tumorAltContexts = tumorAltContexts;
    }

    @NotNull
    public AltContext normal() {
        return normal;
    }

    @NotNull
    public List<AltContext> tumorAltContexts() {
        return tumorAltContexts;
    }

    @NotNull
    public AltContext primaryTumor() {
        return tumorAltContexts.get(0);
    }

}
