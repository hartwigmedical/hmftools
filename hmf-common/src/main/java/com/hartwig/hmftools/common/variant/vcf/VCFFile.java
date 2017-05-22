package com.hartwig.hmftools.common.variant.vcf;

import com.hartwig.hmftools.common.variant.Variant;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public abstract class VCFFile<T extends Variant> {

    @NotNull
    abstract List<String> originalMetaInformationLines();

    @NotNull
    abstract String originalHeaderLine();

    @NotNull
    public abstract List<T> variants();
}
