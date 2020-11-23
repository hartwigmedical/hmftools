package com.hartwig.hmftools.common.clinvar;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ClinvarSummary {

    @NotNull
    String info();

    @NotNull
    ClinvarPathogenicity pathogenicity();
}
