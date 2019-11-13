package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.math.Doubles;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

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

    public double minorAllelePloidy() {
        return Doubles.lessThan(averageActualBAF(), 0.50) ? 0 : Math.max(0, (1 - averageActualBAF()) * averageTumorCopyNumber());
    }

    public double majorAllelePloidy() {
        return averageTumorCopyNumber() - minorAllelePloidy();
    }

    public long length() {
        return end() - start() + 1;
    }

    public abstract long minStart();
    
    public abstract long maxStart();
}