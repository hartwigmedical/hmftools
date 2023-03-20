package com.hartwig.hmftools.datamodel.hla;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LilacAllele {

    @NotNull
    String allele();

    double tumorCopyNumber();

    double somaticMissense();
    double somaticNonsenseOrFrameshift();
    double somaticSplice();
    double somaticSynonymous();
    double somaticInframeIndel();

    // used to create the report
    int refFragments();
    int tumorFragments();
    int rnaFragments();
}
