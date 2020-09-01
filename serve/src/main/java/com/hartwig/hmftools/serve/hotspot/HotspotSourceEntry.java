package com.hartwig.hmftools.serve.hotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface HotspotSourceEntry {

    @NotNull
    String gene();

    @Nullable
    String transcript();

    @NotNull
    String proteinAnnotation();
}
