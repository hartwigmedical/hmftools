package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class GroupBoundary
{
    public final String Chromosome;
    public final int Position;
    public final Orientation Orient;

    public GroupBoundary(final String chromosome, final int position, final Orientation orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orient = orientation;
    }

    public String toString() { return format("%s:%d:%d", Chromosome, Position, Orient.asByte()); }
}
