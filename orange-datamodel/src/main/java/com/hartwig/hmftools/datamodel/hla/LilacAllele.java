package com.hartwig.hmftools.datamodel.hla;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LilacAllele
{
    @NotNull
    String geneClass();

    @NotNull
    String allele();

    @NotNull
    String qcStatus();

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
