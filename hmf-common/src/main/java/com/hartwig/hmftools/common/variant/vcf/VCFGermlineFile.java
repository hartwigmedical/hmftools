package com.hartwig.hmftools.common.variant.vcf;

import com.hartwig.hmftools.common.variant.GermlineVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VCFGermlineFile extends VCFFile<GermlineVariant> {

    @NotNull
    public abstract String refSample();

    @Nullable
    public abstract String tumorSample();
}
