package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._12;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._13;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._14;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._15;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._16;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._20;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._21;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._22;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;

import org.junit.Test;

public class ChromosomeArmCopyNumberTest
{
    @Test
    public void headingsTest()
    {
        String expected = "chromosome\tarm\tmeanCopyNumber\tmedianCopyNumber\tminCopyNumber\tmaxCopyNumber";
        assertEquals(expected, ChromosomeArmCopyNumber.tsvFileHeader());
    }

    @Test
    public void toTsvTest()
    {
        ChromosomeArmCopyNumber o1 = new ChromosomeArmCopyNumber(
                _1, Arm.P, 0.663333333333, 0.600000001, 0.320092, 0.90002000200025);
        assertEquals("1\tP\t0.6633\t0.6000\t0.3201\t0.9000", o1.toTSV());
    }

    @Test
    public void fromTsvTest()
    {
        ChromosomeArmCopyNumber original = new ChromosomeArmCopyNumber(
                _1, Arm.P, 0.6633, 0.6000, 0.3201, 0.9000);
        ChromosomeArmCopyNumber rebuilt = ChromosomeArmCopyNumber.fromTsv(original.toTSV());
        assertEquals(original, rebuilt);
    }

    @Test
    public void isReportable()
    {
        assertTrue(cn(_1, Arm.P).includeInReport());
        assertTrue(cn(_1, Arm.Q).includeInReport());
        assertTrue(cn(_12, Arm.P).includeInReport());
        assertTrue(cn(_12, Arm.Q).includeInReport());
        assertFalse(cn(_13, Arm.P).includeInReport());
        assertTrue(cn(_13, Arm.Q).includeInReport());
        assertFalse(cn(_14, Arm.P).includeInReport());
        assertTrue(cn(_14, Arm.Q).includeInReport());
        assertFalse(cn(_15, Arm.P).includeInReport());
        assertTrue(cn(_15, Arm.Q).includeInReport());
        assertTrue(cn(_16, Arm.P).includeInReport());
        assertTrue(cn(_16, Arm.Q).includeInReport());
        assertTrue(cn(_20, Arm.P).includeInReport());
        assertTrue(cn(_20, Arm.Q).includeInReport());
        assertFalse(cn(_21, Arm.P).includeInReport());
        assertTrue(cn(_21, Arm.Q).includeInReport());
        assertFalse(cn(_22, Arm.P).includeInReport());
        assertTrue(cn(_22, Arm.Q).includeInReport());
        assertTrue(cn(_X, Arm.P).includeInReport());
        assertTrue(cn(_X, Arm.Q).includeInReport());
        assertTrue(cn(_Y, Arm.P).includeInReport());
        assertTrue(cn(_Y, Arm.Q).includeInReport());
    }

    private ChromosomeArmCopyNumber cn(final HumanChromosome chromosome, final Arm arm)
    {
        return new ChromosomeArmCopyNumber(chromosome, arm, 0.0, 0.0, 0.0, 0.0);
    }
}
