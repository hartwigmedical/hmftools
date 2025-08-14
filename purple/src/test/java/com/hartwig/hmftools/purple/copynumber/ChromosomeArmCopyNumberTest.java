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
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ChromosomeArm;

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
        ChromosomeArmCopyNumber o1 = new ChromosomeArmCopyNumber(_1, P_ARM, 0.663333333333, 0.600000001, 0.320092, 0.90002000200025);
        assertEquals("1\tP\t0.6633\t0.6000\t0.3201\t0.9000", o1.toTSV());
    }

    @Test
    public void isReportable()
    {
        assertTrue(cn(_1, P_ARM).includeInReport());
        assertTrue(cn(_1, Q_ARM).includeInReport());
        assertTrue(cn(_12, P_ARM).includeInReport());
        assertTrue(cn(_12, Q_ARM).includeInReport());
        assertFalse(cn(_13, P_ARM).includeInReport());
        assertTrue(cn(_13, Q_ARM).includeInReport());
        assertFalse(cn(_14, P_ARM).includeInReport());
        assertTrue(cn(_14, Q_ARM).includeInReport());
        assertFalse(cn(_15, P_ARM).includeInReport());
        assertTrue(cn(_15, Q_ARM).includeInReport());
        assertTrue(cn(_16, P_ARM).includeInReport());
        assertTrue(cn(_16, Q_ARM).includeInReport());
        assertTrue(cn(_20, P_ARM).includeInReport());
        assertTrue(cn(_20, Q_ARM).includeInReport());
        assertFalse(cn(_21, P_ARM).includeInReport());
        assertTrue(cn(_21, Q_ARM).includeInReport());
        assertFalse(cn(_22, P_ARM).includeInReport());
        assertTrue(cn(_22, Q_ARM).includeInReport());
        assertTrue(cn(_X, P_ARM).includeInReport());
        assertTrue(cn(_X, Q_ARM).includeInReport());
        assertTrue(cn(_Y, P_ARM).includeInReport());
        assertTrue(cn(_Y, Q_ARM).includeInReport());
    }

    private ChromosomeArmCopyNumber cn(HumanChromosome chromosome, ChromosomeArm arm)
    {
        return new ChromosomeArmCopyNumber(chromosome, arm, 0.0, 0.0, 0.0, 0.0);
    }
}
