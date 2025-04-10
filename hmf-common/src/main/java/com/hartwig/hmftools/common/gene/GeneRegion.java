package com.hartwig.hmftools.common.gene;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class GeneRegion extends ChrBaseRegion
{
    public final String GeneName;
    public final int ExonRank;

    public GeneRegion(final String chromosome, final int posStart, final int posEnd, final String geneName, final int exonRank)
    {
        super(chromosome, posStart, posEnd);
        GeneName = geneName;
        ExonRank = exonRank;
    }

    public String toString() { return format("%s exon(%d) region(%s)", GeneName, ExonRank, super.toString()); }
}
