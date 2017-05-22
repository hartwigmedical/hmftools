package com.hartwig.hmftools.common.variant.vcf;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VCFSomaticFile extends VCFFile<SomaticVariant> {

    @NotNull
    public abstract String sample();
}
