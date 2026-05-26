package com.hartwig.hmftools.redux.splice.rescue;

import java.util.Objects;

// (chrom, intronStart, intronEnd) — 1-based inclusive intron coordinates. Used as the key for the
// annotated-junction lookup the JunctionRescueResolver consults. intronStart is the first intronic
// base (= last exonic base + 1) and intronEnd is the last intronic base (= next exonic base - 1).
public class ChrIntron
{
    public final String Chromosome;
    public final int IntronStart;
    public final int IntronEnd;

    public ChrIntron(final String chromosome, final int intronStart, final int intronEnd)
    {
        Chromosome = chromosome;
        IntronStart = intronStart;
        IntronEnd = intronEnd;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o) return true;
        if(!(o instanceof ChrIntron)) return false;
        final ChrIntron other = (ChrIntron) o;
        return IntronStart == other.IntronStart
                && IntronEnd == other.IntronEnd
                && Chromosome.equals(other.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, IntronStart, IntronEnd);
    }

    @Override
    public String toString()
    {
        return Chromosome + ":" + IntronStart + "-" + IntronEnd;
    }
}
