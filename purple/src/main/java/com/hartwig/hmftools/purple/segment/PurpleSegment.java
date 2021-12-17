package com.hartwig.hmftools.purple.segment;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleSegment implements GenomeRegion {

    @NotNull
    @Override
    public abstract String chromosome();

    @Override
    public abstract int start();

    @Override
    public abstract int end();

    public abstract boolean ratioSupport();

    public abstract SegmentSupport support();

    public abstract boolean svCluster();

    public abstract int minStart();

    public abstract int maxStart();
}
