package com.hartwig.hmftools.datamodel.hla;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LilacRecord
{
    @NotNull
    String qc();

    @NotNull
    List<LilacAllele> alleles();
}
