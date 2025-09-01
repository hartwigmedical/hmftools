package com.hartwig.hmftools.common.bam.testutilities;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionDepth extends ChrBaseRegion
{
    public final HumanChromosome mChromosomeIndex;
    private final int depth;

    public RegionDepth(final HumanChromosome chromosome, final int posStart, final int posEnd, final int depth)
    {
        super("chr" + chromosome , posStart, posEnd);
        Preconditions.checkArgument(depth >= 0);
        Preconditions.checkArgument(posStart < posEnd);
        this.depth = depth;
        this.mChromosomeIndex = chromosome;
    }

    public BamRegionWriter singleBlockBamRegionWriter()
    {
        int range = end() - start();
        Preconditions.checkArgument( range % 2 == 0);
        int readLength = range / 2;
        BamRegionWriter bamRegionWriter = new  BamRegionWriter(mChromosomeIndex, start(), end());
        bamRegionWriter.setReadLength(readLength);
        bamRegionWriter.setStepLength(range);
        bamRegionWriter.setDepthAtEachStep(depth);
        return bamRegionWriter;
    }

    public BamRegionWriter length100ReadsBamRegionWriter()
    {
        int range = end() - start();
        Preconditions.checkArgument( range % 200 == 0);

        BamRegionWriter bamRegionWriter = new BamRegionWriter(mChromosomeIndex, start(), end());
        bamRegionWriter.setReadLength(100);
        bamRegionWriter.setStepLength(200);
        bamRegionWriter.setDepthAtEachStep(depth);
        return bamRegionWriter;
    }

    public boolean precedes(RegionDepth other)
    {
        return this.chromosome().equals(other.chromosome()) && this.end() <= other.start();
    }
}
