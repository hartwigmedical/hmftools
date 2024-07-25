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
    // TODO better attribute names
    int numberHotspotMutations();

    int numberHotspotSVs();

    int sumSNPAlleleReadCounts();

    int sumTumorVariantFragmentCountsExclSGL();

    int sumBafCount();
}
