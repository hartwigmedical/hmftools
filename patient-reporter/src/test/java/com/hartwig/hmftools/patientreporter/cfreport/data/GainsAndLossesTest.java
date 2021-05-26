package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.protect.cnchromosome.ChromosomeArmKey;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GainsAndLossesTest {

    @Test
    public void canDetermineCopyNumberPArm() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        String chromosome = "1";
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString(chromosome), ChromosomeArm.P_ARM), 1.123);

        assertEquals("1", GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss(chromosome, "p.12")));
    }

    @Test
    public void canDetermineCopyNumberQArm() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        String chromosome = "4";
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString(chromosome), ChromosomeArm.Q_ARM), 4.51);

        assertEquals("5", GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss(chromosome, "q.12")));
    }

    @Test
    public void crashOnCopyNumberDifferentChromosomes() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.P_ARM), 1.123);

        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("2", "p.12")));
    }

    @Test
    public void crashOnCopyNumberUnknownArms() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.UNKNOWN), 2.34);

        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", "p.12")));
    }

    @Test
    public void crashOnCopyNumberDifferentArms() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM), 2.34);
        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", "p.12")));
    }

    @Test(expected = IllegalArgumentException.class)
    public void crashOnDetermineCopyNumberWeirdArm() {
        Map<ChromosomeArmKey, Double> cnPerChromosome = Maps.newHashMap();
        cnPerChromosome.put(new ChromosomeArmKey(HumanChromosome.fromString("1"), ChromosomeArm.Q_ARM), 1.123);
        GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", ".12"));
    }

    @NotNull
    private static ReportableGainLoss testGainLoss(@NotNull String chromosome, @NotNull String chromosomeBand) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(chromosome)
                .chromosomeBand(chromosomeBand)
                .gene(Strings.EMPTY)
                .copies(0)
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .build();
    }
}