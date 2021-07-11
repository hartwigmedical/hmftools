package com.hartwig.hmftools.purple.recovery;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
interface RecoveredVariant {

    @NotNull
    PurpleCopyNumber copyNumber();

    @NotNull
    PurpleCopyNumber prevCopyNumber();

    @NotNull
    VariantContext context();

    @Nullable
    VariantContext mate();

    @NotNull
    StructuralVariant variant();
}
