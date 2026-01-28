package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._14;
import static com.hartwig.hmftools.common.segmentation.Arm.P;
import static com.hartwig.hmftools.common.segmentation.Arm.Q;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Test;

public class ChrArmLocatorTest
{
    @Test
    public void defaultLocatorTest()
    {
        ChrArmLocator defaultLocator = ChrArmLocator.defaultLocator(RefGenomeVersion.V38);
        assertEquals(new ChrArm(_14, P), defaultLocator.map("chr14", 10_000_000));
        assertEquals(new ChrArm(_14, P), defaultLocator.map("chr14", 17086761));
        assertEquals(new ChrArm(_14, Q), defaultLocator.map("chr14", 17086762));
    }

    @Test
    public void mapWithGenomePositionTest()
    {
        ChrArmLocator defaultLocator = ChrArmLocator.defaultLocator(RefGenomeVersion.V38);
        assertEquals(new ChrArm(_14, P), defaultLocator.map(new GenomePositionImpl("chr14", 10_000_000)));
    }
}
