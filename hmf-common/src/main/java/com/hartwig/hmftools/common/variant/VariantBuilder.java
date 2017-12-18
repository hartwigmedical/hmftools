package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public interface VariantBuilder {
    VariantBuilder type(VariantType type);

    VariantBuilder chromosome(@NotNull final String chromosome);

    VariantBuilder position(final long position);

    VariantBuilder ref(@NotNull final String ref);

    VariantBuilder alt(@NotNull final String alt);

    VariantBuilder filter(@NotNull final String filter);

    SomaticVariant build();
}
