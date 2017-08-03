package com.hartwig.hmftools.common.variant.vcf;

import java.util.List;

import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StrelkaVCFSomaticFile extends VCFFile<StrelkaSomaticVariant> {
    @NotNull
    @Override
    public abstract List<String> originalMetaInformationLines();

    @NotNull
    @Override
    public abstract String originalHeaderLine();

    @NotNull
    @Override
    public abstract List<StrelkaSomaticVariant> variants();
}
