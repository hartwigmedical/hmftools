package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.hartwig.hmftools.common.position.GenomeInterval;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Segment of DNA
 */
public class BgSegment implements GenomeInterval {
    private PurpleCopyNumber cn;
    public BgSegment(PurpleCopyNumber cn) {
        this.cn = cn;
    }
    public static BgSegment createUnplacedSegment() {
        return new BgSegment(ImmutablePurpleCopyNumber.builder()
                .segmentEndSupport(SegmentSupport.NONE)
                .segmentStartSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .chromosome("unplaced")
                .start(1)
                .end(1)
                .averageTumorCopyNumber(0)
                .averageActualBAF(0)
                .averageObservedBAF(0)
                .build());
    }
    public double ploidy() {
        return cn.averageTumorCopyNumber();
    }
    @Nullable
    @Override
    public Integer startOffset() {
        return 0;
    }

    @Nullable
    @Override
    public Integer endOffset() {
        return (int)(cn.end() - cn.start());
    }

    @NotNull
    @Override
    public String chromosome() {
        return cn.chromosome();
    }

    @Override
    public long position() {
        return cn.start();
    }
    private static boolean canMerge(BgSegment left, BgSegment right) {
        return left.cn.end() == right.cn.start() - 1 && left.cn.chromosome().equals(right.cn.chromosome());
    }
    public static BgSegment merge(BgSegment left, BgSegment right) {
        if (!canMerge(left, right)) {
            throw new IllegalArgumentException("Cannot merge segments");
        }
        if (left.cn.segmentEndSupport() == SegmentSupport.CENTROMERE ||  right.cn.segmentStartSupport() == SegmentSupport.CENTROMERE) {
            throw new IllegalArgumentException("Merging across centromere");
        }
        BgSegment merged = new BgSegment(merge(left.cn, right.cn));
        return merged;
    }
    private static PurpleCopyNumber merge(PurpleCopyNumber left, PurpleCopyNumber right) {
        long length = left.length() + right.length();
        int bafCount = left.bafCount() + right.bafCount();
        ImmutablePurpleCopyNumber merged = ImmutablePurpleCopyNumber.builder()
                .chromosome(left.chromosome())
                .start(left.start())
                .end(right.end())
                .segmentStartSupport(left.segmentStartSupport())
                .segmentEndSupport(right.segmentEndSupport())
                .averageTumorCopyNumber(
                        left.averageTumorCopyNumber() * left.length() / (double) length +
                                right.averageTumorCopyNumber() * right.length() / (double) length)
                .bafCount(bafCount)
                .averageActualBAF(
                        left.averageActualBAF() * left.bafCount() / (double) bafCount +
                                right.averageActualBAF() * right.bafCount() / (double) bafCount)
                .depthWindowCount(left.depthWindowCount() + right.depthWindowCount())
                .method(left.method() == right.method() ? left.method() : CopyNumberMethod.UNKNOWN)
                .averageObservedBAF(left.averageObservedBAF() * left.bafCount() / (double) bafCount +
                        right.averageObservedBAF() * right.bafCount() / (double) bafCount)
                .build();
        return merged;
    }
}
