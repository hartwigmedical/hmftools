package com.hartwig.hmftools.datamodel.purple;

import org.jetbrains.annotations.NotNull;

// TODO: remove this redundant interface after refactoring `ReportableVariant` in OncoAct.
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
