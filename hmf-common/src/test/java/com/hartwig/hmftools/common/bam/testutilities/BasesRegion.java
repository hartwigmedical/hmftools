package com.hartwig.hmftools.common.bam.testutilities;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class BasesRegion extends ChrBaseRegion
{
    public final byte[] mBases;
    public final HumanChromosome mChromosome;

    public BasesRegion(final HumanChromosome chromosome, int start, int end, byte[] bases)
    {
        super(V38.versionedChromosome(chromosome), start, end);
        mBases = bases;
        mChromosome = chromosome;
    }
}
