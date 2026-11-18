package com.hartwig.hmftools.amber.blacklist;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._5;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;

import org.junit.Before;
import org.junit.Test;

public class AmberBlacklistStatisticsTest
{
    private AmberBlacklistStatistics statistics;

    @Before
    public void setup()
    {
        statistics = new AmberBlacklistStatistics();
    }

    @Test
    public void returnedPointsAreSortedByPosition()
    {
        record(_3, 2000, 0.2, 20);
        record(_1, 2000, 0.2, 20);
        record(_2, 2000, 0.2, 20);
        record(_2, 1000, 0.2, 20);
        record(_1, 1000, 0.2, 20);
        record(_3, 1000, 0.2, 20);
        List<AmberBlacklistPoint> result = statistics.findSuspiciousPoints(20);
        assertEquals(6, result.size());
        checkPosition(_1, 1000, result.get(0));
        checkPosition(_1, 2000, result.get(1));
        checkPosition(_2, 1000, result.get(2));
        checkPosition(_2, 2000, result.get(3));
        checkPosition(_3, 1000, result.get(4));
        checkPosition(_3, 2000, result.get(5));
    }

    @Test
    public void positionsNotRepresentedByAtLeastAFifthOfSamplesAreNotReturned()
    {
        record(_3, 3000, 0.1, 18);
        record(_3, 2000, 0.3, 20);
        record(_3, 1000, 0.2, 21);
        record(_2, 2000, 0.2, 25);
        record(_1, 2000, 0.2, 30);

        List<AmberBlacklistPoint> result = statistics.findSuspiciousPoints(100);
        assertEquals(3, result.size());
        checkPosition(_1, 2000, result.get(0));
        checkPosition(_2, 2000, result.get(1));
        checkPosition(_3, 1000, result.get(2));
    }

    @Test
    public void positionsWithAverageVafInNormalBandAreNotReturned()
    {
        record(_2, 4000, 0.39, 10);
        record(_2, 4000, 0.41, 11);

        record(_3, 4000, 0.39, 11);
        record(_3, 4000, 0.41, 10);

        record(_4, 4000, 0.59, 10);
        record(_4, 4000, 0.61, 11);

        record(_5, 4000, 0.59, 11);
        record(_5, 4000, 0.61, 10);

        List<AmberBlacklistPoint> result = statistics.findSuspiciousPoints(25);
        assertEquals(2, result.size());
        checkPosition(_3, 4000, result.get(0));
        checkPosition(_4, 4000, result.get(1));
    }

    @Test
    public void countsTest()
    {
        record(_3, 3000, 0.1, 100);
        record(_3, 4000, 0.2, 123);
        List<AmberBlacklistPoint> result = statistics.findSuspiciousPoints(125);
        assertEquals(2, result.size());
        assertEquals(100, result.get(0).count());
        assertEquals(123, result.get(1).count());
    }

    @Test
    public void meanVafTest()
    {
        record(_3, 3000, 0.1, 100);
        record(_3, 3000, 0.2, 100);
        record(_3, 4000, 0.2, 120);
        record(_3, 4000, 0.6, 30);
        List<AmberBlacklistPoint> result = statistics.findSuspiciousPoints(300);
        assertEquals(2, result.size());
        assertEquals(0.15, result.get(0).meanVaf(), 0.0001);
        assertEquals(0.28, result.get(1).meanVaf(), 0.0001);
    }

    private void checkPosition(HumanChromosome chromosome, int position, AmberBlacklistPoint statistic)
    {
        assertEquals(V38.versionedChromosome(chromosome), statistic.chromosome());
        assertEquals(position, statistic.position());
    }

    private void record(HumanChromosome chromosome, int position, double vaf, int count)
    {
        for(int i = 0; i < count; i++)
        {
            statistics.record(gp(chromosome, position), vaf);
        }
    }

    private GenomePosition gp(HumanChromosome chromosome, int position)
    {
        return new GenomePositionImpl(V38.versionedChromosome(chromosome), position);
    }
}
