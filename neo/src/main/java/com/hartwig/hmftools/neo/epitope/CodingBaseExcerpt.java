package com.hartwig.hmftools.neo.epitope;

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
}
