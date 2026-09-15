package com.hartwig.hmftools.purple.reseg;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class Supersegment
{
    // A contiguous run of 'bothNone' segments (unsupported by an SV on either side, diploid), potentially
    // bridging over isolated germline events
    public final int Id;
    public final List<ObservedRegion> BothNoneMembers;
    public final List<ObservedRegion> SkippableMembers;

    public Supersegment(final int id, final List<ObservedRegion> bothNoneMembers, final List<ObservedRegion> skippableMembers)
    {
        Id = id;
        BothNoneMembers = bothNoneMembers;
        SkippableMembers = skippableMembers;
    }

    public String chromosome() { return BothNoneMembers.get(0).chromosome(); }
}
