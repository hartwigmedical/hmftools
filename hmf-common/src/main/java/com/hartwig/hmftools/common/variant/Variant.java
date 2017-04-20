package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public interface Variant {

    @NotNull
    VariantType type();

    @NotNull
    String filter();
}
