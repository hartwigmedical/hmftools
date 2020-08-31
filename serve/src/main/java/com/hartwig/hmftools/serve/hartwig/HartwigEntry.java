package com.hartwig.hmftools.serve.hartwig;

import com.hartwig.hmftools.serve.hotspot.HotspotSourceEntry;

import org.jetbrains.annotations.NotNull;

public interface HartwigEntry extends HotspotSourceEntry {

    @NotNull
    String chromosome();

    long position();

    @NotNull
    String ref();

    @NotNull
    String alt();
}
