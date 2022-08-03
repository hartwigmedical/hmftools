package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

public class GroupBoundary
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

    public GroupBoundary(final String chromosome, final int position, final byte orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
    }

    public String toString() { return format("%s:%d:%d", Chromosome, Position, Orientation); }
}
