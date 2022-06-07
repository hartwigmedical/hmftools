package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.cnchromosome.ImmutableCnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.interpretation.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.common.utils.DataUtil;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GainsAndLossesTest {

    @Test
    public void canDetermineCopyNumberPArm() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        String chromosome = "1";

        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString(chromosome))
                .chromosomeArm(ChromosomeArm.P_ARM)
                .copyNumber(1.123).build());

        assertEquals("1", GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss(chromosome, "p.12")));
    }

    @Test
    public void canDetermineCopyNumberQArm() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        String chromosome = "4";
        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString(chromosome))
                .chromosomeArm(ChromosomeArm.Q_ARM)
                .copyNumber(4.51).build());

        assertEquals("5", GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss(chromosome, "q.12")));
    }

    @Test
    public void crashOnCopyNumberDifferentChromosomes() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString("1"))
                .chromosomeArm(ChromosomeArm.P_ARM)
                .copyNumber(1.123).build());

        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("2", "p.12")));
    }

    @Test
    public void crashOnCopyNumberUnknownArms() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString("1"))
                .chromosomeArm(ChromosomeArm.UNKNOWN)
                .copyNumber(2.34).build());
        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", "p.12")));
    }

    @Test
    public void crashOnCopyNumberDifferentArms() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString("1"))
                .chromosomeArm(ChromosomeArm.Q_ARM)
                .copyNumber(2.34).build());
        assertEquals(DataUtil.NA_STRING, GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", "p.12")));
    }

    @Test(expected = IllegalArgumentException.class)
    public void crashOnDetermineCopyNumberWeirdArm() {
        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        cnPerChromosome.add(ImmutableCnPerChromosomeArmData.builder()
                .chromosome(HumanChromosome.fromString("1"))
                .chromosomeArm(ChromosomeArm.Q_ARM)
                .copyNumber(1.123).build());
        GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, testGainLoss("1", ".12"));
    }

    @NotNull
    private static ReportableGainLoss testGainLoss(@NotNull String chromosome, @NotNull String chromosomeBand) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(chromosome)
                .chromosomeBand(chromosomeBand)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false)
                .minCopies(0)
                .maxCopies(0)
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .build();
    }
}