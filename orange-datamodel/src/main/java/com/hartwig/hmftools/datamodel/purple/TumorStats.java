package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TumorStats
{
    double maxDiploidProportion();

    int hotspotMutationCount();

    int hotspotStructuralVariantCount();

    int smallVariantAlleleReadCount();

    int structuralVariantTumorFragmentCount();

    int bafCount();
}
