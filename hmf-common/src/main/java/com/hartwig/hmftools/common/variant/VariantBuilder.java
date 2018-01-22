package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public interface VariantBuilder {
    VariantBuilder type(@NotNull VariantType type);

    VariantBuilder chromosome(@NotNull String chromosome);

    VariantBuilder position(long position);

    VariantBuilder ref(@NotNull String ref);

    VariantBuilder alt(@NotNull String alt);

    VariantBuilder filter(@NotNull String filter);

    SomaticVariant build();
}
