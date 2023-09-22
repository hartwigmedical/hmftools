package com.hartwig.hmftools.common.genome.refgenome;

import static java.lang.String.format;

import java.util.Comparator;

public class CoordMapping
{
    public final String Chromosome;
    public final int SourceStart;
    public final int SourceEnd;
    public final int DestStart;
    public final int DestEnd;
    public final boolean Reverse;

    public CoordMapping(
            final String chromosome, final int sourceStart, final int sourceEnd, final int destStart, final int destEnd,
            final boolean reverse)
    {
        Chromosome = chromosome;
        SourceStart = sourceStart;
        SourceEnd = sourceEnd;
        DestStart = destStart;
        DestEnd = destEnd;
        Reverse = reverse;
    }

    public int convertPosition(int position)
    {
        int posDiff = position - SourceStart;

        if(!Reverse)
            return DestStart + posDiff;
        else
            return DestEnd - posDiff;
    }

    public int reversePosition(int position)
    {
        int posDiff = !Reverse ? position - DestStart : DestEnd - position;
        return posDiff + SourceStart;
    }

    public String toString()
    {
        return format("%s:%d-%d -> %d-%d %s", Chromosome, SourceStart, SourceEnd, DestStart, DestEnd, Reverse ? "reversed" : "");
    }

    public static class CoordComparatorOnDestination implements Comparator<CoordMapping>
    {
        public int compare(final CoordMapping first, final CoordMapping second)
        {
            if(first.DestStart == second.DestStart)
                return 0;

            return first.DestStart < second.DestStart ? -1 : 1;
        }
    }

}
