package com.hartwig.hmftools.sage.phase;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

public class MnvCandidate implements VariantHotspot {
    private final VariantHotspot mnv;
    private final int lps;

    public MnvCandidate(final VariantHotspot mnv, final int lps) {
        this.mnv = mnv;
        this.lps = lps;
    }

    public VariantHotspot mnv() {
        return mnv;
    }

    public int lps() {
        return lps;
    }

    @NotNull
    @Override
    public String ref() {
        return mnv.ref();
    }

    @NotNull
    @Override
    public String alt() {
        return mnv.alt();
    }

    @NotNull
    @Override
    public String chromosome() {
        return mnv.chromosome();
    }

    @Override
    public long position() {
        return mnv.position();
    }
}

