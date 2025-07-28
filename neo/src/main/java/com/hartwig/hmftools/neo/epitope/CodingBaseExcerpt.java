package com.hartwig.hmftools.neo.epitope;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import htsjdk.samtools.Cigar;

public class CodingBaseExcerpt
{
    public final String Bases;
    public final int[] Positions;
    public final Cigar CigarRef;

    public CodingBaseExcerpt(final String bases, final int upPos, final int downPos, final Cigar cigarRef)
    {
        Bases = bases;
        Positions = new int[] { upPos, downPos};
        CigarRef = cigarRef;
    }

    public String toString()
    {
        return format("baseLen(%d) pos(%d - %d) cigar(%s)",
            Bases.length(), Positions[SE_START], Positions[SE_END], CigarRef.toString());
    }
}
