package com.hartwig.hmftools.panelbuilder;

import java.util.Comparator;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public record NamedRegion(
        ChrBaseRegion region,
        String name
) implements Comparable<NamedRegion>
{
    @Override
    public int compareTo(@NotNull final NamedRegion other)
    {
        return Comparator.comparing(NamedRegion::region).thenComparing(NamedRegion::name).compare(this, other);
    }
}
