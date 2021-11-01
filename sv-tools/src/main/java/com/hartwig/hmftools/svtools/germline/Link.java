package com.hartwig.hmftools.svtools.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;

public class Link
{
    public final String VcfId;
    public final String OtherVcfId;
    public final int MinDistance;
    public final int MaxDistance;

    public Link(final String vcfId, final String otherVcfId, final int minDistance, final int maxDistance)
    {
        VcfId = vcfId;
        OtherVcfId = otherVcfId;
        MinDistance = minDistance;
        MaxDistance = maxDistance;
    }

    public static Link from(final SvData sv)
    {
        int duplicationLength = sv.type() == DUP ? sv.length() + 1 : 0;
        int distance = duplicationLength + sv.insertSequenceLength();

        return new Link(sv.contextStart().getID(), sv.contextEnd().getID(), distance, distance);
    }

    public static Link from(final Breakend first, final Breakend second)
    {
        int minDistance = abs(first.maxPosition() - second.minPosition());
        int maxDistance = abs(first.minPosition() - second.maxPosition());

        return new Link(first.VcfId, second.VcfId, minDistance, maxDistance);
    }

    public String toString() { return String.format("%s - %s distance(%d - %d)", VcfId, OtherVcfId, MinDistance, MaxDistance); }
}
