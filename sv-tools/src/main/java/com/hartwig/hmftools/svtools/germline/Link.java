package com.hartwig.hmftools.svtools.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

public class Link
{
    public final String Id;

    private final Breakend[] mBreakends;
    public final int MinDistance;
    public final int MaxDistance;

    public Link(final String id, final Breakend first, final Breakend second, final int minDistance, final int maxDistance)
    {
        Id = id;

        mBreakends = new Breakend[] { first, second };
        MinDistance = minDistance;
        MaxDistance = maxDistance;
    }

    public static Link from(final SvData sv)
    {
        int distance = sv.duplicationLength() + sv.insertSequenceLength();
        return new Link("SV", sv.breakendStart(), sv.breakendEnd(), distance, distance);
    }

    public static Link from(final String linkId, final Breakend first, final Breakend second)
    {
        int minDistance = abs(first.maxPosition() - second.minPosition());
        int maxDistance = abs(first.minPosition() - second.maxPosition());

        return new Link(linkId, first, second, minDistance, maxDistance);
    }

    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public Breakend otherBreakend(final Breakend breakend) { return breakendStart() == breakend ?  breakendEnd() : breakendStart(); }

    public String toString() { return String.format("%s<%s>%s distance(%d - %d)",
            breakendStart().VcfId, breakendEnd().VcfId, MinDistance, MaxDistance); }
}
