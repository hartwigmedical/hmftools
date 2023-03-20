package com.hartwig.hmftools.datamodel.hla;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public interface LilacRecord {
    @NotNull
    String qc();

    @NotNull
    List<LilacAllele> alleles();
}
