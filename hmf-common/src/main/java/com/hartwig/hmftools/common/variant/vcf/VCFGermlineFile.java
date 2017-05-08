package com.hartwig.hmftools.common.variant.vcf;

import com.hartwig.hmftools.common.variant.GermlineVariant;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VCFGermlineFile extends VCFFile<GermlineVariant> {

    @NotNull
    public abstract String refSample();

    @NotNull
    public abstract String tumorSample();
}
