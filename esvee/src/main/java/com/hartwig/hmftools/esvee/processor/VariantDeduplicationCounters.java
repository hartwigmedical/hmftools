package com.hartwig.hmftools.esvee.processor;

import com.hartwig.hmftools.esvee.util.Counter;

public class VariantDeduplicationCounters
{
    public final Counter VariantsRemoved = new Counter("Variants Removed", false);
    public final Counter OverallTimeNanos = new Counter("Variant Deduplication Time", true);

    public final Counter MatchedGroupWallTimeNanos = new Counter("Matched Descriptor Wall Time", true);
    public final Counter MatchedGroupItemTimeNanos = new Counter("Matched Descriptor Item Time", true);
    public final Counter MatchedGroupRemoved = new Counter("Identical Variants Removed", false);

    public final Counter AssemblyDedupeWallTimeNanos = new Counter("Assembly Dedupe Wall Time", true);
    public final Counter AssemblyDedupeItemTimeNanos = new Counter("Assembly Dedupe Item Time", true);
    public final Counter AssembliesRemoved = new Counter("Assemblies Removed", false);

    public final Counter SinglesDedupeWallTimeNanos = new Counter("SGL Dedupe Wall Time", true);
    public final Counter SinglesDedupeItemTimeNanos = new Counter("SGL Dedupe Item Time", true);
    public final Counter SinglesDeduped = new Counter("SGL Variants Removed", false);

    public final Counter NearbyDedupeWallTimeNanos = new Counter("Nearby Dedupe Wall Time", true);
    public final Counter NearbyDedupeItemTimeNanos = new Counter("Nearby Dedupe Item Time", true);
    public final Counter NearbyDeduped = new Counter("Variants Merged /w Neighbours", false);
}
