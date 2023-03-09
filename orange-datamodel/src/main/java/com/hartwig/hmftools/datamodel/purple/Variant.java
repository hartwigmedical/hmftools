package com.hartwig.hmftools.datamodel.purple;

import org.jetbrains.annotations.NotNull;

public interface Variant {

    @NotNull
    PurpleVariantType type();

    @NotNull
    String gene();

    @NotNull
    String chromosome();

    int position();

    @NotNull
    String ref();

    @NotNull
    String alt();
}
