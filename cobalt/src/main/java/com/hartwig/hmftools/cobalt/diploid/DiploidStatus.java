package com.hartwig.hmftools.cobalt.diploid;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class DiploidStatus extends ChrBaseRegion
{
    public final boolean isDiploid;
    public DiploidStatus(final String chromosome, final int posStart, final int posEnd, final boolean isDiploid)
    {
        super(chromosome, posStart, posEnd);
        this.isDiploid = isDiploid;
    }

    @Override
    public String toString()
    {
        return String.format("%s:%d-%d %b", Chromosome, start(), end(), isDiploid);
    }
}
