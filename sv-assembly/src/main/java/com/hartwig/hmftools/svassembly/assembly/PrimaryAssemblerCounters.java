package com.hartwig.hmftools.svassembly.assembly;

import com.hartwig.hmftools.svassembly.util.Counter;

public class PrimaryAssemblerCounters extends Counters<PrimaryAssemblerCounters>
{
    private static final String CATEGORY = PrimaryAssemblerCounters.class.getSimpleName();

    public final Counter ProcessTimeNanos = new Counter(CATEGORY, "Process Time", true);
    public final Counter InitialReadTimeNanos = new Counter(CATEGORY, "Initial Read Time", true);
    public final Counter DiagramGenerationTimeNanos = new Counter(CATEGORY, "Diagram Generation Time", true);
    public final Counter GraphSimplificationTimeNanos = new Counter(CATEGORY, "Graph Simplification Time", true);
    public final Counter JunctionConstructionTimeNanos = new Counter(CATEGORY, "Junction Construction", true);
    public final Counter JunctionExtensionTimeNanos = new Counter(CATEGORY, "Junction Extension", true);
    public final Counter AnchorConstructionTimeNanos = new Counter(CATEGORY, "Junction Anchoring", true);

    public final Counter ReadsCrossingJunction = new Counter(CATEGORY, "# Reads crossing junction");
    public final Counter ReadsPassingRawQualityThreshold = new Counter(CATEGORY, "# Reads >baseq threshold");
    public final Counter ReadsPassingJunctionQualityThreshold = new Counter(CATEGORY, "# Reads >baseq threshold at junction");
    public final Counter ReadsSoftClippedAtJunction = new Counter(CATEGORY, "# Reads soft-clipped at junction");
    public final Counter HasAcceptableMapQ = new Counter(CATEGORY, "# Reads with acceptable MapQ");
    public final Counter WellMapped = new Counter(CATEGORY, "# Reads not badly mapped");

    public final Counter FlattenedInitial = new Counter(CATEGORY, "Flattened Initial");
    public final Counter InitialAssemblies = new Counter(CATEGORY, "Initial Assemblies");
    public final Counter DedupedInitialAssemblies = new Counter(CATEGORY, "Deduped Initial Assemblies");
    public final Counter FlattenedAnchors = new Counter(CATEGORY, "Flattened Anchors");
    public final Counter AnchoredAssemblies = new Counter(CATEGORY, "Anchored Assemblies");
    public final Counter DedupedAnchoredAssemblies = new Counter(CATEGORY, "Deduped Anchored Assemblies");
}
