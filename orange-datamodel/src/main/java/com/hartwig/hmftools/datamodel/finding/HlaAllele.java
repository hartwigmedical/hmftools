package com.hartwig.hmftools.datamodel.finding;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface HlaAllele extends Finding
{
    @NotNull
    String allele();

    double tumorCopyNumber();

    @Nullable
    Integer refFragments();

    int tumorFragments();

    @Nullable
    Integer rnaFragments();

    double somaticMissense();

    double somaticNonsenseOrFrameshift();

    double somaticSplice();

    double somaticSynonymous();

    double somaticInframeIndel();
}
