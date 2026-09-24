package com.hartwig.hmftools.purple.reseg;

import java.util.List;

import com.hartwig.hmftools.purple.region.ObservedRegion;

public class SupersegmentResult
{
    public final List<Supersegment> Supersegments;
    public final List<ObservedRegion> ExcludedSegments;

    public SupersegmentResult(final List<Supersegment> supersegments, final List<ObservedRegion> excludedSegments)
    {
        Supersegments = supersegments;
        ExcludedSegments = excludedSegments;
    }
}
