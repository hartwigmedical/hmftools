package com.hartwig.hmftools.common.variant.vcf;

import java.util.List;

import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public abstract class VCFFile<T extends Variant> {

    @NotNull
    abstract List<String> originalMetaInformationLines();

    @NotNull
    abstract String originalHeaderLine();

    @NotNull
    public abstract List<T> variants();
}
