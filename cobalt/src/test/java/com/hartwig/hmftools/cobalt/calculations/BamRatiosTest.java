package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.when;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.function.Function;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class BamRatiosTest extends CalculationsTestBase
{
    private final Random RandomGenerator = new Random();
    private BamRatios mBamRatios;
    private Map<Chromosome, List<Integer>> OriginalPositions;
    private Map<Chromosome, List<Double>> OriginalRatios;
    private Map<Chromosome, List<Double>> OriginalDepths;
    private Map<Chromosome, List<Double>> OriginalGCs;

    @Before
    public void setup()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for (int i=0; i<100; i++)
        {
            int position = i * 1000 + 1;
            ratios.put(_1, randomRatio(_1, position));
            ratios.put(_2, randomRatio(_2, position));
            ratios.put(_3, randomRatio(_3, position));
        }
        mBamRatios = new BamRatios(ratios);
        OriginalPositions = extractPositions();
        OriginalRatios = extractRatios();
        OriginalGCs = extractGCs();
        OriginalDepths = extractDepths();
    }

    @Test
    public void normaliseTest()
    {
        class RecordingNormaliser implements ResultsNormaliser
        {
            final Set<BamRatio> ratiosRecorded = new HashSet<>();
            final Set<BamRatio> ratiosNormalised = new HashSet<>();
            boolean DataCollectionFinishedCalled = false;

            @Override
            public void recordValue(final BamRatio bamRatio)
            {
                ratiosRecorded.add(bamRatio);
            }

            @Override
            public void dataCollectionFinished()
            {
                // This should be called after collection but before normalisation.
                assertEquals(ratiosRecorded.size(), mBamRatios.Ratios.size());
                assertTrue(ratiosNormalised.isEmpty());
                DataCollectionFinishedCalled = true;
            }

            @Override
            public void normalise(final BamRatio bamRatio)
            {
                ratiosNormalised.add(bamRatio);
            }
        }
        RecordingNormaliser normaliser = new RecordingNormaliser();
        mBamRatios.normalise(normaliser);

        assertEquals(mBamRatios.Ratios.values().size(), normaliser.ratiosRecorded.size());
        assertEquals(mBamRatios.Ratios.values().size(), normaliser.ratiosNormalised.size());
        for (BamRatio ratio : mBamRatios.Ratios.values())
        {
            Assert.assertTrue(normaliser.ratiosRecorded.contains(ratio));
            Assert.assertTrue(normaliser.ratiosNormalised.contains(ratio));
        }
        assertTrue(normaliser.DataCollectionFinishedCalled);
    }

    @Test
    public void consolidationIdentical()
    {
        mBamRatios.consolidate(new NoOpConsolidator());
        assertEquals(OriginalPositions, extractPositions());
        assertEquals(OriginalRatios, extractRatios());
        assertEquals(OriginalGCs, extractGCs());
        assertEquals(OriginalDepths, extractDepths());
    }

    @Test
    public void consolidateBy10()
    {
        ListMultimap<Chromosome, BamRatio> consolidatedRatios = ArrayListMultimap.create();
        for (int i=1; i<=10; i++)
        {
            int position = i * 10_000 + 1;
            consolidatedRatios.put(_1, randomRatio(_1, position));
            consolidatedRatios.put(_2, randomRatio(_2, position));
            consolidatedRatios.put(_3, randomRatio(_3, position));
        }
        ResultsConsolidator consolidator = Mockito.mock(ResultsConsolidator.class);
        when(consolidator.consolidate(mBamRatios.Ratios)).thenReturn(consolidatedRatios);

        mBamRatios.consolidate(consolidator);

        // The positions, depths and gc values should be unchanged.
        assertEquals(OriginalPositions, extractPositions());
        assertEquals(OriginalGCs, extractGCs());
        assertEquals(OriginalDepths, extractDepths());

        // The ratios are -1.0 if the position does not match that of a consolidated ratio,
        // and the value of the corresponding consolidated ratio if it does.
        for (int i=1; i<100; i++)
        {
            int position = i * 1000 + 1;
            if (i % 10 == 0)
            {
                int c = (i / 10) - 1;
                assertEquals(position,mBamRatios.Ratios.get(_1).get(i).position());
                assertEquals(position, consolidatedRatios.get(_1).get(c).position());
                assertEquals(consolidatedRatios.get(_1).get(c).ratio(),mBamRatios.Ratios.get(_1).get(i).ratio(),  0.0001);
                assertEquals(consolidatedRatios.get(_2).get(c).ratio(),mBamRatios.Ratios.get(_2).get(i).ratio(),  0.0001);
                assertEquals(consolidatedRatios.get(_3).get(c).ratio(),mBamRatios.Ratios.get(_3).get(i).ratio(),  0.0001);
            }
            else
            {
                assertEquals(-1.0,mBamRatios.Ratios.get(_1).get(i).ratio(),  0.0001);
                assertEquals(-1.0,mBamRatios.Ratios.get(_2).get(i).ratio(),  0.0001);
                assertEquals(-1.0,mBamRatios.Ratios.get(_3).get(i).ratio(),  0.0001);
            }
        }
    }

    @Test
    public void consolidateBy30()
    {
        ListMultimap<Chromosome, BamRatio> consolidatedRatios = ArrayListMultimap.create();
        for (int i=1; i<=3; i++)
        {
            int position = i * 30_000 + 1;
            consolidatedRatios.put(_1, randomRatio(_1, position));
            consolidatedRatios.put(_2, randomRatio(_2, position));
            consolidatedRatios.put(_3, randomRatio(_3, position));
        }
        ResultsConsolidator consolidator = Mockito.mock(ResultsConsolidator.class);
        when(consolidator.consolidate(mBamRatios.Ratios)).thenReturn(consolidatedRatios);

        mBamRatios.consolidate(consolidator);

        // The positions, depths and gc values should be unchanged.
        assertEquals(OriginalPositions, extractPositions());
        assertEquals(OriginalGCs, extractGCs());
        assertEquals(OriginalDepths, extractDepths());

        // The ratios are -1.0 if the position does not match that of a consolidated ratio,
        // and the value of the corresponding consolidated ratio if it does.
        for (int i=1; i<100; i++)
        {
            int position = i * 1000 + 1;
            if (i % 30 == 0)
            {
                int c = (i / 30) - 1;
                assertEquals(position,mBamRatios.Ratios.get(_1).get(i).position());
                assertEquals(position, consolidatedRatios.get(_1).get(c).position());
                assertEquals(consolidatedRatios.get(_1).get(c).ratio(),mBamRatios.Ratios.get(_1).get(i).ratio(),  0.0001);
                assertEquals(consolidatedRatios.get(_2).get(c).ratio(),mBamRatios.Ratios.get(_2).get(i).ratio(),  0.0001);
                assertEquals(consolidatedRatios.get(_3).get(c).ratio(),mBamRatios.Ratios.get(_3).get(i).ratio(),  0.0001);
            }
            else
            {
                assertEquals(-1.0,mBamRatios.Ratios.get(_1).get(i).ratio(),  0.0001);
                assertEquals(-1.0,mBamRatios.Ratios.get(_2).get(i).ratio(),  0.0001);
                assertEquals(-1.0,mBamRatios.Ratios.get(_3).get(i).ratio(),  0.0001);
            }
        }
    }

    @Test
    public void consolidateEmpty()
    {
        mBamRatios = new BamRatios(ArrayListMultimap.create());
        ListMultimap<Chromosome, BamRatio> consolidatedRatios = ArrayListMultimap.create();
        for (int i=1; i<=5; i++)
        {
            int position = i * 10_000 + 1;
            consolidatedRatios.put(_1, randomRatio(_1, position));
        }
        ResultsConsolidator consolidator = Mockito.mock(ResultsConsolidator.class);
        when(consolidator.consolidate(mBamRatios.Ratios)).thenReturn(consolidatedRatios);

        mBamRatios.consolidate(consolidator);
        assertEquals(0, mBamRatios.Ratios.size());
    }

    BamRatio randomRatio(Chromosome chromosome, int position)
    {
        double ratio = RandomGenerator.nextDouble();
        double depth = RandomGenerator.nextDouble() * 10.0;
        double gc = (RandomGenerator.nextDouble() * 40.0 + 25)/100.1;
        return new BamRatio(chromosome, position, depth, ratio, gc);
    }

    private Map<Chromosome, List<Integer>> extractPositions()
    {
        return  extractValues(BamRatio::position);
    }

    private Map<Chromosome, List<Double>> extractRatios()
    {
        return extractValues(BamRatio::ratio);
    }

    private Map<Chromosome, List<Double>> extractGCs()
    {
        return extractValues(BamRatio::gcContent);
    }

    private Map<Chromosome, List<Double>> extractDepths()
    {
        return extractValues(BamRatio::readDepth);
    }

    private <T extends Number> Map<Chromosome, List<T>> extractValues(Function<BamRatio, T> f)
    {
        Map<Chromosome, List<T>> result = new HashMap<>();
        mBamRatios.Ratios.keySet().forEach(chromosome -> {
            List<T> chrPositions = mBamRatios.Ratios.get(chromosome).stream().map(f).toList();
            result.put(chromosome, chrPositions);
        });
        return result;
    }
}
