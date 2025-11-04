package com.hartwig.hmftools.cobalt.normalisers;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;

import java.util.List;

import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;

import org.junit.Assert;
import org.junit.Test;

public class DoNothingNormaliserTest extends CalculationsTestBase
{
    @Test
    public void recordTest()
    {
        BamRatio br1 = br(_1, 1001, 10.0, 0.45, true);
        BamRatio br2 = br(_1, 2001, 20.0, 0.46, true);
        BamRatio br3 = br(_3, 2001, 30.0, 0.46, true);
        List<BamRatio> bamRatios = of(br1, br2, br3);
        DoNothingNormaliser normaliser = new DoNothingNormaliser();
        bamRatios.forEach(normaliser::recordValue);
        Assert.assertEquals(10.0, br1.ratio(), 0.001);
        Assert.assertEquals(20.0, br2.ratio(), 0.001);
        Assert.assertEquals(30.0, br3.ratio(), 0.001);

        bamRatios.forEach(normaliser::normalise);
        Assert.assertEquals(10.0, br1.ratio(), 0.001);
        Assert.assertEquals(20.0, br2.ratio(), 0.001);
        Assert.assertEquals(30.0, br3.ratio(), 0.001);
    }
}
