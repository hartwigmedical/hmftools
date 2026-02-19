package com.hartwig.hmftools.common.genome.position;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._14;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import org.junit.Assert;
import org.junit.Test;

public class GenomePositionImplTest
{
    @Test
    public void fromChromosomeConstructorTest()
    {
        GenomePositionImpl gpi = new GenomePositionImpl(_14, V38, 89);
        assertEquals(89, gpi.position());
        assertEquals(V38.versionedChromosome(_14), gpi.chromosome());
    }

    @Test
    public void copyConstructorTest()
    {
        GenomePositionImpl original = new GenomePositionImpl(_14, V38, 89);
        GenomePositionImpl copy = new GenomePositionImpl(original);

        assertEquals(original.position(), copy.position());
        assertEquals(original.chromosome(), copy.chromosome());
        Assert.assertNotSame(original, copy);
        assertEquals(original, copy);
    }
}
