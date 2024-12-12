package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;

public class UnmapRegion extends ChrBaseRegion
{
    public final HighDepthRegion Region;
    public final boolean IsStandardChromosome; // distinction recorded but not currently used for any specific logic

    public static final UnmapRegion UNMAPPED_READS = new UnmapRegion(NO_CHROMOSOME_NAME, NO_POSITION, NO_POSITION, null);

    public UnmapRegion(final String chromosome, final int posStart, final int posEnd, final HighDepthRegion region)
    {
        super(chromosome, posStart, posEnd);
        IsStandardChromosome = HumanChromosome.contains(chromosome);
        Region = region;
    }

    public String toString()
    {
        return format("%s depth(%d)", super.toString(), Region.maxDepth());
    }
}
