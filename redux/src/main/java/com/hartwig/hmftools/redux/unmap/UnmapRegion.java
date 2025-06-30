package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;

public class UnmapRegion extends ChrBaseRegion
{
    public final List<UnmappingRegion> Regions;
    public final boolean IsStandardChromosome; // distinction recorded but not currently used for any specific logic

    public static final UnmapRegion UNMAPPED_READS = new UnmapRegion(NO_CHROMOSOME_NAME, NO_POSITION, NO_POSITION, null);

    public UnmapRegion(final String chromosome, final int posStart, final int posEnd, final UnmappingRegion region)
    {
        super(chromosome, posStart, posEnd);
        IsStandardChromosome = HumanChromosome.contains(chromosome);
        Regions = Lists.newArrayList(region);
    }

    public String toString()
    {
        return format("%s depth(%d)", super.toString(), Regions.stream().mapToInt(x -> x.maxDepth()).max().orElse(0));
    }
}
