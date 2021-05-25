package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.protect.cnchromosome.ChromosomeArmKey;

import org.junit.Test;

public class GainsAndLossesTest {

    @Test
    public void testCanDetermineCopyNumberParm() {

        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM), 1.123);
        String chromosome = "1";
        String chromosomeBand = "p.12";

        assertEquals("1", GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand));
    }

    @Test
    public void testCanDetermineCopyNumberQarm() {

        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("4"), ChromosomeArm.Q_ARM), 4.51);
        String chromosome = "4";
        String chromosomeBand = "q.12";

        assertEquals("5", GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand));
    }

    @Test(expected = NullPointerException.class)
    public void crashTestCanDetermineCopyNumberDiffferentChromosomes() {

        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM), 1.123);
        String chromosome = "2";
        String chromosomeBand = "p.12";

        GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand);
    }

//    @Test(expected = NullPointerException.class)
//    public void crashTestCanDetermineCopyNumberUnknownArms() {
//
//        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
//        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.UNKNOWN), 2.34);
//        String chromosome = "1";
//        String chromosomeBand = "p.12";
//
//        GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand);
//    }
//
//    @Test(expected = NullPointerException.class)
//    public void crashTestCanDetermineCopyNumberDifferentArms() {
//
//        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
//        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM), 2.34);
//        String chromosome = "1";
//        String chromosomeBand = "p.12";
//
//        GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand);
//    }
//
//    @Test(expected = NullPointerException.class)
//    public void crashTestCanDetermineCopyNumber() {
//        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
//        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM), 1.123);
//        String chromosome = "1";
//        String chromosomeBand = ".12";
//
//       GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand);
//    }
}