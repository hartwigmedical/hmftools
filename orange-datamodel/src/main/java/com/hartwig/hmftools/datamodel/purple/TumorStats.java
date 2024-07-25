package com.hartwig.hmftools.datamodel.purple;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorStats {
    // TODO better attribute names
    public abstract int numberHotspotMutations();

    public abstract int numberHotspotSVs();

    public abstract int sumSNPAlleleReadCounts();

    public abstract int sumTumorVariantFragmentCountsExclSGL();

    public abstract int sumBafCount();
}
