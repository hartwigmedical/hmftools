package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.genome.region.GenomeRegion;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements GenomeRegion {

    public abstract int bafCount();

    public abstract double averageActualBAF();

    public abstract double averageObservedBAF();

    public abstract double averageTumorCopyNumber();

    public abstract int depthWindowCount();

    public abstract SegmentSupport segmentStartSupport();

    public abstract SegmentSupport segmentEndSupport();

    public abstract CopyNumberMethod method();

    public abstract double gcContent();

    public abstract int minStart();

    public abstract int maxStart();
}