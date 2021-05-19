package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReporter;
import com.hartwig.hmftools.protect.cnchromosome.ChromosomeArmKey;

import org.junit.Test;

public class GainsAndLossesTest {

    @Test
    public void testCanDetermineCopyNumber() {

        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM), 1.123);
        String chromosome = "1";
        String chromosomeBand = "p.12";

        assertEquals("1", GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand));

    }

    @Test(expected = IllegalStateException.class)
    public void crashTestCanDetermineCopyNumber() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM), 1.123);
        String chromosome = "1";
        String chromosomeBand = ".12";

       GainsAndLosses.copyChromosomeArm(cnPerChromosome, chromosome, chromosomeBand);

    }

}