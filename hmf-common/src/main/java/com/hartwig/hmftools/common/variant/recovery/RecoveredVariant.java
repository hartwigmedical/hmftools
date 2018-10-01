package com.hartwig.hmftools.common.variant.recovery;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RecoveredVariant extends GenomeRegion {

    double baf();

    double copyNumber();

    int depthWindowCount();

    double prevBaf();

    double prevCopyNumber();

    long prevLength();

    int prevDepthWindowCount();

    SegmentSupport previous();

    SegmentSupport support();

    SegmentSupport next();


    @Nullable
    String variant();

    @Nullable
    String mate();

    @Nullable
    Double qual();

    @Nullable
    String filter();

    @Nullable
    Integer orientation();

    @Nullable
    String alt();

}
