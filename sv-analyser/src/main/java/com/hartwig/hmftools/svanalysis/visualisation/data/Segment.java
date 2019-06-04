package com.hartwig.hmftools.svanalysis.visualisation.data;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.svanalysis.visualisation.circos.SegmentTerminal;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Segment implements GenomeRegion {

    public abstract String sampleId();

    public abstract int clusterId();

    public abstract int chainId();

    public abstract int track();

    public abstract int traverseCount();

    public abstract SegmentTerminal startTerminal();

    public abstract SegmentTerminal endTerminal();
}
