package com.hartwig.hmftools.datamodel.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class LilacSummaryData {
    @NotNull
    public abstract String qc();

    @NotNull
    public abstract List<LilacAllele> alleles();

    public int somaticVariantCount() {
        return (int) alleles().stream().mapToDouble(LilacAllele::somaticVariantCount).sum();
    }
}
