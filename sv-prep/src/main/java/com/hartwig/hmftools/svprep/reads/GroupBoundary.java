package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import com.hartwig.hmftools.common.sv.Direction;

public class GroupBoundary
{
    public final String Chromosome;
    public final int Position;
    public final Direction Orientation;

    public GroupBoundary(final String chromosome, final int position, final Direction orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
    }

    public String toString() { return format("%s:%d:%d", Chromosome, Position, Orientation.Step); }
}
