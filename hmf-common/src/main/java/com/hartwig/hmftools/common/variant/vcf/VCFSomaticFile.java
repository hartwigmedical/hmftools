package com.hartwig.hmftools.common.variant.vcf;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VCFSomaticFile extends VCFFile<SomaticVariant> {

    @NotNull
    public abstract String sample();
}
