package com.hartwig.hmftools.serve.hotspot;

import org.jetbrains.annotations.NotNull;

public interface HotspotSourceEntry {

    @NotNull
    String gene();

    @NotNull
    String transcript();

    @NotNull
    String proteinAnnotation();
}
