package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public interface Variant {

    @NotNull
    String chromosome();

    long position();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    VariantType type();

    @NotNull
    String filter();
}
