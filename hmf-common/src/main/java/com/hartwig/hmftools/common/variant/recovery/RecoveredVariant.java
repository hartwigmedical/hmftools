package com.hartwig.hmftools.common.variant.recovery;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RecoveredVariant extends GenomeRegion {

    // Copy Number
    double baf();

    double copyNumber();

    int depthWindowCount();

    double gcContent();

    double prevBaf();

    double prevCopyNumber();

    long prevLength();

    int prevDepthWindowCount();

    double prevGCContent();

    SegmentSupport previous();

    SegmentSupport support();

    SegmentSupport next();

    double nextGCContent();

    long minStart();

    long maxStart();

    @Nullable
    String variant();

    @Nullable
    Double qual();

    @Nullable
    String filter();

    @Nullable
    Integer orientation();

    @Nullable
    String mate();

    @Nullable
    Integer mateOrientation();

    @Nullable
    Long mateMinStart();

    @Nullable
    Long mateMaxStart();

    @Nullable
    SegmentSupport mateSupport();

    @Nullable
    String alt();

}
