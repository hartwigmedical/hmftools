package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.cobalt.count.ReadDepth;

import org.junit.Test;

public class BamRatiosMeanCalculatorTest
{

    @Test
    public void mean()
    {
        ReadDepth rd1 = new ReadDepth("1", 1001, 80, 0.49);
        BamRatio ratio1 = new BamRatio(_1, rd1, true);
        ratio1.normaliseForGc(100.0);

        ReadDepth rd2 = new ReadDepth("1", 1001, 80, 0.46);
        BamRatio ratio2 = new BamRatio(_1, rd2, true);
        ratio2.normaliseForGc(80);

        ReadDepth rd3 = new ReadDepth("2", 1001, 90, 0.59);
        BamRatio ratio3 = new BamRatio(_2, rd3, true);
        ratio3.normaliseForGc(100.0);

        ReadDepth rd4 = new ReadDepth("2", 1001, 800, 0.49);
        BamRatio ratio4 = new BamRatio(_2, rd4, true);
        ratio4.normaliseForGc(-1.0);

        BamRatiosMeanCalculator calculator = new BamRatiosMeanCalculator();
        calculator.recordValue(ratio1);
        calculator.recordValue(ratio2);
        calculator.recordValue(ratio3);
        calculator.recordValue(ratio4); // not included

        assertEquals(0.9, calculator.mean(), 0.001);
    }

    @Test
    public void skipZeroValues()
    {
        ReadDepth rd1 = new ReadDepth("1", 1001, 0.0, 0.49);
        BamRatio ratio1 = new BamRatio(_1, rd1, true);
        ratio1.normaliseForGc(100.0);

        ReadDepth rd2 = new ReadDepth("1", 1001, 80, 0.46);
        BamRatio ratio2 = new BamRatio(_1, rd2, true);
        ratio2.normaliseForGc(80);

        BamRatiosMeanCalculator calculator = new BamRatiosMeanCalculator();
        calculator.recordValue(ratio1);
        calculator.recordValue(ratio2);
        assertEquals(1.0, calculator.mean(), 0.001);
    }

    @Test
    public void excludeAllosomes()
    {
        ReadDepth rd1 = new ReadDepth("1", 1001, 200.0, 0.49);
        BamRatio ratio1 = new BamRatio(_1, rd1, true);
        ratio1.normaliseForGc(100.0);

        ReadDepth rd2 = new ReadDepth("1", 1001, 80, 0.46);
        BamRatio ratio2 = new BamRatio(_1, rd2, true);
        ratio2.normaliseForGc(80);

        ReadDepth rd3 = new ReadDepth("X", 1001, 50, 0.46);
        BamRatio ratio3 = new BamRatio(_X, rd3, true);
        ratio3.normaliseForGc(10);

        ReadDepth rd4 = new ReadDepth("Y", 1001, 60, 0.46);
        BamRatio ratio4 = new BamRatio(_X, rd4, true);
        ratio4.normaliseForGc(10);

        BamRatiosMeanCalculator calculator = new BamRatiosMeanCalculator();
        calculator.recordValue(ratio1);
        calculator.recordValue(ratio2);
        calculator.recordValue(ratio3);
        calculator.recordValue(ratio4);
        assertEquals(1.5, calculator.mean(), 0.001);
    }
}
