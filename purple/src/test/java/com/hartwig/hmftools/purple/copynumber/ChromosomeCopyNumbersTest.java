package com.hartwig.hmftools.purple.copynumber;

import static java.util.List.of;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createCopyNumber;
import static com.hartwig.hmftools.common.purple.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.TELOMERE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

public class ChromosomeCopyNumbersTest
{
    private static final double EPSILON = 0.000000001;

    @Test
    public void onePurpleCopyNumber()
    {
        final double copyNumber = 0.55;
        PurpleCopyNumber pcn1 = pcn(_1, 10_000_000, 100_000_000, copyNumber, TELOMERE);
        ChromosomeCopyNumbers copyNumbers = new ChromosomeCopyNumbers(of(pcn1));
        List<ChromosomeArmCopyNumber> results = copyNumbers.data();
        assertEquals(1, results.size());
        ChromosomeArmCopyNumber result = results.get(0);
        assertEquals(_1, result.chromosome());
        assertEquals(P_ARM, result.arm());
        assertEquals(copyNumber, result.meanCopyNumber(), EPSILON);
        assertEquals(copyNumber, result.medianCopyNumber(), EPSILON);
        assertEquals(copyNumber, result.maxCopyNumber(), EPSILON);
        assertEquals(copyNumber, result.minCopyNumber(), EPSILON);
    }

    @Test
    public void oneArm()
    {
        PurpleCopyNumber pcn1 = pcn(_1, 1, 10_000_000, 0.5, TELOMERE);
        PurpleCopyNumber pcn2 = pcn(_1, 10_000_001, 20_000_000, 0.5, NONE);
        PurpleCopyNumber pcn3 = pcn(_1, 20_000_001, 30_000_000, 0.5, NONE);
        PurpleCopyNumber pcn4 = pcn(_1, 30_000_001, 40_000_000, 0.7, NONE);
        PurpleCopyNumber pcn5 = pcn(_1, 40_000_001, 50_000_000, 0.8, NONE);

        ChromosomeCopyNumbers copyNumbers = new ChromosomeCopyNumbers(of(pcn1, pcn2, pcn3, pcn4, pcn5));
        List<ChromosomeArmCopyNumber> results = copyNumbers.data();
        assertEquals(1, results.size());
        ChromosomeArmCopyNumber result = results.get(0);
        assertEquals(_1, result.chromosome());
        assertEquals(P_ARM, result.arm());
        assertEquals(0.6, result.meanCopyNumber(), EPSILON);
        assertEquals(0.5, result.medianCopyNumber(), EPSILON);
        assertEquals(0.8, result.maxCopyNumber(), EPSILON);
        assertEquals(0.5, result.minCopyNumber(), EPSILON);
    }

    @Test
    public void meanIsWeightedBySegmentLength()
    {
        PurpleCopyNumber pcn1 = pcn(_1, 1, 100, 0.4, NONE); // 100 * 0.4 = 40
        PurpleCopyNumber pcn2 = pcn(_1, 101, 400, 0.5, NONE); // 300 * 0.5 = 150
        PurpleCopyNumber pcn3 = pcn(_1, 401, 500, 0.6, NONE); // 100 * 0.6 = 60
        PurpleCopyNumber pcn4 = pcn(_1, 501, 1000, 0.7, NONE); // 500 * 0.7 = 350
        PurpleCopyNumber pcn5 = pcn(_1, 1001, 2000, 0.8, NONE); // 1000 * 0.8 = 800

        ChromosomeCopyNumbers copyNumbers = new ChromosomeCopyNumbers(of(pcn1, pcn2, pcn3, pcn4, pcn5));
        List<ChromosomeArmCopyNumber> results = copyNumbers.data();
        ChromosomeArmCopyNumber result = results.get(0);
        // 40 + 150 + 60 + 350 + 800 = 1400, total length = 2000, so weighted average = 0.7
        assertEquals(0.7, result.meanCopyNumber(), EPSILON);
        // There are 2000 positions, so the median is the 1000th
        assertEquals(0.8, result.medianCopyNumber(), EPSILON);
    }

    @Test
    public void medianIsWeightedBySegmentLength()
    {
        PurpleCopyNumber pcn1 = pcn(_1, 1, 100, 0.4, NONE); // 100 @ 0.4
        PurpleCopyNumber pcn2 = pcn(_1, 101, 400, 0.5, NONE); // 300 @ 0.5
        PurpleCopyNumber pcn3 = pcn(_1, 401, 600, 0.4, NONE); // 200 @ 0.4
        PurpleCopyNumber pcn4 = pcn(_1, 601, 1000, 0.7, NONE); // 400 @ 0.7
        PurpleCopyNumber pcn5 = pcn(_1, 1001, 1200, 0.3, NONE); // 200 @ 0.3
        PurpleCopyNumber pcn6 = pcn(_1, 1201, 1600, 0.7, NONE); // 400 @ 0.7
        PurpleCopyNumber pcn7 = pcn(_1, 1601, 1700, 0.6, NONE); // 100 @ 0.6
        PurpleCopyNumber pcn8 = pcn(_1, 1701, 1900, 0.8, NONE); // 200 @ 0.8

        ChromosomeCopyNumbers copyNumbers = new ChromosomeCopyNumbers(of(pcn1, pcn2, pcn3, pcn4, pcn5, pcn6, pcn7, pcn8));
        List<ChromosomeArmCopyNumber> results = copyNumbers.data();
        ChromosomeArmCopyNumber result = results.get(0);
        // 200 100 200 300 100 400 400 200, total is 1700, median value is 950th
        assertEquals(0.7, result.medianCopyNumber(), EPSILON);
    }

    @Test
    public void multipleArms()
    {
        PurpleCopyNumber pcn1p1 = pcn(_1, 1, 50_000_000, 0.5, TELOMERE);
        PurpleCopyNumber pcn1p2 = pcn(_1, 50_000_001, 100_000_000, 0.7, NONE);
        PurpleCopyNumber pcn1q1 = pcn(_1, 160_000_001, 200_000_000, 0.4, CENTROMERE);
        PurpleCopyNumber pcn1q2 = pcn(_1, 200_000_001, 240_000_000, 0.6, NONE);
        PurpleCopyNumber pcn2p1 = pcn(_2, 1, 40_000_000, 0.55, TELOMERE);
        PurpleCopyNumber pcn2p2 = pcn(_2, 40_000_001, 80_000_000, 0.75, NONE);
        PurpleCopyNumber pcn2q1 = pcn(_2, 100_000_001, 150_000_000, 0.45, CENTROMERE);
        PurpleCopyNumber pcn2q2 = pcn(_2, 150_000_001, 200_000_000, 0.65, NONE);
        ChromosomeCopyNumbers copyNumbers =
                new ChromosomeCopyNumbers(of(pcn1p1, pcn1p2, pcn1q1, pcn1q2, pcn2p1, pcn2p2, pcn2q1, pcn2q2));
        List<ChromosomeArmCopyNumber> results = copyNumbers.data();
        assertEquals(4, results.size());
        ChromosomeArmCopyNumber p1 = results.get(0);
        assertEquals(_1, p1.chromosome());
        assertEquals(P_ARM, p1.arm());
        assertEquals(0.6, p1.meanCopyNumber(), EPSILON);
        assertEquals(0.7, p1.medianCopyNumber(), EPSILON);
        assertEquals(0.7, p1.maxCopyNumber(), EPSILON);
        assertEquals(0.5, p1.minCopyNumber(), EPSILON);

        ChromosomeArmCopyNumber q1 = results.get(1);
        assertEquals(_1, q1.chromosome());
        assertEquals(Q_ARM, q1.arm());
        assertEquals(0.5, q1.meanCopyNumber(), EPSILON);
        assertEquals(0.6, q1.medianCopyNumber(), EPSILON);
        assertEquals(0.6, q1.maxCopyNumber(), EPSILON);
        assertEquals(0.4, q1.minCopyNumber(), EPSILON);

        ChromosomeArmCopyNumber p2 = results.get(2);
        assertEquals(_2, p2.chromosome());
        assertEquals(P_ARM, p2.arm());
        assertEquals(0.65, p2.meanCopyNumber(), EPSILON);
        assertEquals(0.75, p2.medianCopyNumber(), EPSILON);
        assertEquals(0.75, p2.maxCopyNumber(), EPSILON);
        assertEquals(0.55, p2.minCopyNumber(), EPSILON);

        ChromosomeArmCopyNumber q2 = results.get(3);
        assertEquals(_2, q2.chromosome());
        assertEquals(Q_ARM, q2.arm());
        assertEquals(0.55, q2.meanCopyNumber(), EPSILON);
        assertEquals(0.65, q2.medianCopyNumber(), EPSILON);
        assertEquals(0.65, q2.maxCopyNumber(), EPSILON);
        assertEquals(0.45, q2.minCopyNumber(), EPSILON);
    }

    @Test
    public void excludeNonReportableArms()
    {
        List<PurpleCopyNumber> copyNumbers = new ArrayList<>();
        for (HumanChromosome humanChromosome : HumanChromosome.values())
        {
            copyNumbers.add(pcn(humanChromosome, 1, 100, 0.5, NONE)); // P
            copyNumbers.add(pcn(humanChromosome, 1, 100, 0.5, CENTROMERE)); // Q
        }
        ChromosomeCopyNumbers ccn = new ChromosomeCopyNumbers(copyNumbers);
        List<ChromosomeArmCopyNumber> results = ccn.data();
        results.forEach(copyNumber -> assertTrue(copyNumber.includeInReport()));
    }

    private PurpleCopyNumber pcn(HumanChromosome chromosome, int start, int end, double copyNumber, SegmentSupport startSupport)
    {
        return createCopyNumber(chromosome.name().substring(1), start, end, copyNumber, startSupport, SegmentSupport.NONE).build();
    }
}
