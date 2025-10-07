package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class UnityNormaliserTest extends CalculationsTestBase
{
    @Test
    public void recordTest()
    {
        BamRatio br1 = br(_1, 1001, 10.0, 0.45, true);
        BamRatio br2 = br(_1, 2001, 20.0, 0.46, false);
        BamRatio br3 = br(_3, 2001, 30.0, 0.46, true);
        List<BamRatio> bamRatios = of(br1, br2, br3);
        UnityNormaliser normaliser = new UnityNormaliser();
        bamRatios.forEach(normaliser::recordValue);
        Assert.assertEquals(10.0, br1.ratio(), 0.001);
        Assert.assertEquals(-1.0, br2.ratio(), 0.001);
        Assert.assertEquals(30.0, br3.ratio(), 0.001);

        bamRatios.forEach(normaliser::normalise);
        Assert.assertEquals(10.0/20.0, br1.ratio(), 0.001);
        Assert.assertEquals(-1.0, br2.ratio(), 0.001);
        Assert.assertEquals(30.0/20.0, br3.ratio(), 0.001);
    }

    @Test
    public void allosomesDoNotContributeToNormalisationValue()
    {
        BamRatio br1 = br(_1, 1001, 1.0, 0.45, true);
        BamRatio br2 = br(_X, 2001, 2.0, 0.46, true);
        BamRatio br3 = br(_Y, 2001, 3.0, 0.46, true);
        List<BamRatio> bamRatios = of(br1, br2, br3);
        UnityNormaliser normaliser = new UnityNormaliser();
        bamRatios.forEach(normaliser::recordValue);
        Assert.assertEquals(1.0, br1.ratio(), 0.001);
        Assert.assertEquals(2.0, br2.ratio(), 0.001);
        Assert.assertEquals(3.0, br3.ratio(), 0.001);

        bamRatios.forEach(normaliser::normalise);
        Assert.assertEquals(1.0, br1.ratio(), 0.001);
        Assert.assertEquals(2.0, br2.ratio(), 0.001);
        Assert.assertEquals(3.0, br3.ratio(), 0.001);
    }

    @Test
    public void onlyIncludePositiveValues()
    {
        BamRatio br1 = br(_1, 1001, 10.0, 0.45, true);
        BamRatio br2 = br(_1, 2001, 0.0, 0.46, true);
        BamRatio br3 = br(_3, 2001, 50.0, 0.46, true);
        BamRatio br4 = br(_3, 2001, 10.0, 0.46, true);
        br4.applyEnrichment(-1.0); // sets the ratio to be -1.0
        List<BamRatio> bamRatios = of(br1, br2, br3, br4);
        UnityNormaliser normaliser = new UnityNormaliser();
        bamRatios.forEach(normaliser::recordValue);
        Assert.assertEquals(10.0, br1.ratio(), 0.001);
        Assert.assertEquals(0.0, br2.ratio(), 0.001);
        Assert.assertEquals(50.0, br3.ratio(), 0.001);
        Assert.assertEquals(-1.0, br4.ratio(), 0.001);

        bamRatios.forEach(normaliser::normalise);
        Assert.assertEquals(10.0/30.0, br1.ratio(), 0.001);
        Assert.assertEquals(-1.0, br2.ratio(), 0.001);
        Assert.assertEquals(50.0/30.0, br3.ratio(), 0.001);
        Assert.assertEquals(-1.0, br4.ratio(), 0.001);
    }
}
