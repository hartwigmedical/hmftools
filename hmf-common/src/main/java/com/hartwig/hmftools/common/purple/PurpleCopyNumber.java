package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;

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

    public double minorAlleleCopyNumber() {
        return Doubles.lessThan(averageActualBAF(), 0.50) ? 0 : Math.max(0, (1 - averageActualBAF()) * averageTumorCopyNumber());
    }

    public double majorAlleleCopyNumber() {
        return averageTumorCopyNumber() - minorAlleleCopyNumber();
    }

    public int length() { return end() - start() + 1; }

    public abstract int minStart();

    public abstract int maxStart();

    public boolean svSupport() {
        return segmentStartSupport().isSV() || segmentEndSupport().isSV();
    }
}