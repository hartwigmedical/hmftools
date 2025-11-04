package com.hartwig.hmftools.cobalt.normalisers;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Before;
import org.junit.Test;

public class DiploidNormaliserTest extends CalculationsTestBase
{
    private Chromosome CurrentChromosome;
    private int Position = 1;
    private double CurrentGC;
    private DiploidNormaliser normaliser;
    private ListMultimap<Chromosome, BamRatio> ChromosomeToBamRatio;

    @Before
    public void setup()
    {
        CurrentGC = 0.45;
        Position = -1;
        CurrentChromosome = null;
        ChromosomeToBamRatio = ArrayListMultimap.create();
        normaliser = new DiploidNormaliser(5, 7, V38);
    }

    @Test
    public void medianRatiosTest()
    {
        CurrentChromosome = _1;
        Position = 1;
        Pair<Double, Integer> chr1Stats = addRatios(-1, 22, 23, 24, -1, 34, 27, 25, -1, 30, 27, 27, 25, 23, 21, 29, 23, 25, 31, 80, 32);

        CurrentChromosome = _4;
        Position = 1;
        Pair<Double, Integer> chr4Stats = addRatios(19, 19, 20, 22, 23, 19, 17, 20, 27, 27, 23, 33, 21, 29, 23, 25, 31, 23, 23);

        CurrentChromosome = _2;
        Position = 1;
        Pair<Double, Integer> chr2Stats = addRatios(-1, 32, 33, 34, 36, 27, 25, 30, 27, 37, 35, 33, 31, 29, 23, 25, 31, 80, 23);

        normaliser.dataCollectionFinished();

        List<MedianRatio> medianRatios = normaliser.medianRatios();
        assertEquals(3, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(chr1Stats.getLeft(), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals((long) chr1Stats.getRight(), medianRatios.get(0).Count);
        assertEquals(V38.versionedChromosome(_2), medianRatios.get(1).Chromosome);
        assertEquals(chr2Stats.getLeft(), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals((long) chr2Stats.getRight(), medianRatios.get(1).Count);
        assertEquals(V38.versionedChromosome(_4), medianRatios.get(2).Chromosome);
        assertEquals(chr4Stats.getLeft(), medianRatios.get(2).MedianRatio, 0.001);
        assertEquals((long) chr4Stats.getRight(), medianRatios.get(2).Count);
    }

    @Test
    public void normaliseTest()
    {
        normaliser = new DiploidNormaliser(3, 3, V38);

        CurrentChromosome = _1;
        Position = 1;
        Pair<Double, Integer> chr1Stats = addRatios(20, 20, 20, 25, 25, 25);

        CurrentChromosome = _4;
        Position = 1;
        Pair<Double, Integer> chr4Stats = addRatios(20, 20, 20, 20, 30, 30, 30);

        CurrentChromosome = _X;
        Position = 1;
        Pair<Double, Integer> chrXStats = addRatios(-1, 30, 30,50, 30,  30, 40, 40);

        CurrentChromosome = _Y;
        Position = 1;
        Pair<Double, Integer> chrYStats = addRatios(10, 10, 15, -1.0, 10, 15);

        normaliser.dataCollectionFinished();
        normalise();

        List<BamRatio> chr1RatiosNormalised = ChromosomeToBamRatio.get(_1);
        assertEquals(6, chr1RatiosNormalised.size());
        int index = 0;
        assertEquals( 1.0, chr1RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr1RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 20.0/22.5, chr1RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 25.0/22.5, chr1RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr1RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr1RatiosNormalised.get(index).getDiploidAdjustedRatio(), 0.001);

        List<BamRatio> chr4RatiosNormalised = ChromosomeToBamRatio.get(_4);
        assertEquals(7, chr4RatiosNormalised.size());
        index = 0;
        assertEquals( 1.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 30.0/25.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr4RatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 1.0, chr4RatiosNormalised.get(index).getDiploidAdjustedRatio(), 0.001);

        List<BamRatio> chrXRatiosNormalised = ChromosomeToBamRatio.get(_X);
        assertEquals(8, chrXRatiosNormalised.size());
        index = 0;
        assertEquals( -1.0, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5 * 50.0/30.0, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5 * 30.0/35.0, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5, chrXRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.5 * 40.0/35.0, chrXRatiosNormalised.get(index).getDiploidAdjustedRatio(), 0.001);

        // ChrY does not get normalised.
        List<BamRatio> chrYRatiosNormalised = ChromosomeToBamRatio.get(_Y);
        assertEquals(6, chrYRatiosNormalised.size());
        index = 0;
        assertEquals( 0.1, chrYRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.1, chrYRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.15, chrYRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( -1.0, chrYRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.1, chrYRatiosNormalised.get(index++).getDiploidAdjustedRatio(), 0.001);
        assertEquals( 0.15, chrYRatiosNormalised.get(index).getDiploidAdjustedRatio(), 0.001);

        List<MedianRatio> medianRatios = normaliser.medianRatios();
        assertEquals(4, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(chr1Stats.getLeft(), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals((long) chr1Stats.getRight(), medianRatios.get(0).Count);

        assertEquals(V38.versionedChromosome(_4), medianRatios.get(1).Chromosome);
        assertEquals(chr4Stats.getLeft(), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals((long) chr4Stats.getRight(), medianRatios.get(1).Count);

        assertEquals(V38.versionedChromosome(_X), medianRatios.get(2).Chromosome);
        assertEquals(chrXStats.getLeft(), medianRatios.get(2).MedianRatio, 0.001);
        assertEquals((long) chrXStats.getRight(), medianRatios.get(2).Count);

        assertEquals(V38.versionedChromosome(_Y), medianRatios.get(3).Chromosome);
        assertEquals(chrYStats.getLeft(), medianRatios.get(3).MedianRatio, 0.001);
        assertEquals((long) chrYStats.getRight(), medianRatios.get(3).Count);
    }

    private void normalise()
    {
        ChromosomeToBamRatio.keySet().forEach(chromosome ->
        {
            List<BamRatio> chrBamRatios = ChromosomeToBamRatio.get(chromosome);
            chrBamRatios.forEach(chrBamRatio -> normaliser.normalise(chrBamRatio));
        });
    }

    private Pair<Double, Integer> addRatios(Number... depths)
    {
        List<Double> doubles = new ArrayList<>();
        for(Number value : depths)
        {
            BamRatio bamRatio = br(CurrentChromosome, Position, value.doubleValue(), CurrentGC, true);
            if(value.doubleValue() > 0)
            {
                bamRatio.normaliseByMean(100.0);
            }
            doubles.add(value.doubleValue());
            ChromosomeToBamRatio.put(CurrentChromosome, bamRatio);
            normaliser.recordValue(bamRatio);
            Position += 1000;
        }
        final List<Double> positiveValues = doubles.stream().filter(d -> d > 0).collect(Collectors.toList());
        double median = Doubles.median(positiveValues) / 100.0;
        int count = positiveValues.size();

        return new ImmutablePair<>(median, count);
    }
}
