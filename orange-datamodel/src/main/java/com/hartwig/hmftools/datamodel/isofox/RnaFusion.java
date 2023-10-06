package com.hartwig.hmftools.datamodel.isofox;

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
    default String display()
    {
        return String.format("%s_%s", geneStart(), geneEnd());
    }

    @NotNull
    String geneStart();

    @NotNull
    String geneEnd();

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
