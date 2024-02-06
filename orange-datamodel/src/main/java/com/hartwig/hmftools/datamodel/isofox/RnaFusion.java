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
        return String.format("%s::%s", geneStart() != null ? geneStart() : "", geneEnd() != null ? geneEnd() : "");
    }

    @Nullable
    String geneStart();

    @Nullable
    String geneEnd();

    @NotNull
    String chromosomeStart();

    @NotNull
    String chromosomeEnd();

    int positionStart();

    int positionEnd();

    @NotNull
    String junctionTypeStart();

    @NotNull
    String junctionTypeEnd();

    @NotNull
    StructuralVariantType svType();

    int splitFragments();

    int realignedFrags();

    int discordantFrags();

    int depthStart();

    int depthEnd();

    int cohortFrequency();
}
