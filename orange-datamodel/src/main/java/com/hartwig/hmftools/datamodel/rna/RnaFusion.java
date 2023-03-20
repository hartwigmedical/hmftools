package com.hartwig.hmftools.datamodel.rna;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RnaFusion
{

    @NotNull
    String name();

    @NotNull
    String chromosomeUp();

    @NotNull
    String chromosomeDown();

    int positionUp();

    int positionDown();

    @NotNull
    String junctionTypeUp();

    @NotNull
    String junctionTypeDown();

    @NotNull
    StructuralVariantType svType();

    int splitFragments();

    int realignedFrags();

    int discordantFrags();

    int depthUp();

    int depthDown();

    int cohortFrequency();
}
