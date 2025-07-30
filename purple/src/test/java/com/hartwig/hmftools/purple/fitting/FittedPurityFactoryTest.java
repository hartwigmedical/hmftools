package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.cobalt.CobaltTestUtils.male;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.FittingRegion;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.region.ObservedRegionFile;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class FittedPurityFactoryTest extends FittingTestBase
{
    private static final double EPSILON = 1e-10;
    private ExecutorService executorService;
    private RegionFitCalculator mRegionFitCalculator;
    private CobaltChromosomes mCobaltChromosomes;
    private String chr1 = HumanChromosome._1.toString();
    private FittedPurityFactory factory;
    private List<FittedPurity> results;
    private FittedPurity first, second, third;

    @Before
    public void setup()
    {
        super.setup();
        executorService = Executors.newFixedThreadPool(1);
        mCobaltChromosomes = male();
        mRegionFitCalculator = new RegionFitCalculator(mCobaltChromosomes, mConfig.Fitting, 10);
    }

    @After
    public void cleanup()
    {
        executorService.shutdown();
    }

    @Test
    public void testPloidyRange()
    {
        List<Double> oneToThree = FittedPurityFactory.ploidyRange(1, 3);
        assertEquals(101, oneToThree.size());
        assertEquals(1d, oneToThree.get(0), EPSILON);
        assertEquals(3d, oneToThree.get(100), EPSILON);

        List<Double> threeToFive = FittedPurityFactory.ploidyRange(3, 5);
        assertEquals(41, threeToFive.size());
        assertEquals(3d, threeToFive.get(0), EPSILON);
        assertEquals(5d, threeToFive.get(40), EPSILON);

        List<Double> fiveOnwards = FittedPurityFactory.ploidyRange(5, 8);
        assertEquals(31, fiveOnwards.size());
        assertEquals(5d, fiveOnwards.get(0), EPSILON);
        assertEquals(8d, fiveOnwards.get(30), EPSILON);

        List<Double> all = FittedPurityFactory.ploidyRange(1, 8);
        assertEquals(171, all.size());
        assertEquals(1d, all.get(0), EPSILON);
        assertEquals(3d, all.get(100), EPSILON);
        assertEquals(5d, all.get(140), EPSILON);
        assertEquals(8d, all.get(170), EPSILON);
    }

    @Test
    public void testFixedPloidyRange()
    {
        double fixedPloidy = 0.3456;
        List<Double> fixed = FittedPurityFactory.ploidyRange(fixedPloidy, fixedPloidy);
        assertEquals(1, fixed.size());
        assertEquals(fixedPloidy, fixed.get(0), EPSILON);
    }

    @Test
    public void singleObservedRegionNoStructuralVariants() throws Exception
    {
        List<FittingRegion> regions = Lists.newArrayList();
        regions.add(observedRegion(1, 1000, 100, 0.5, 1.0, 0.0));
        createFactoryAndRun(regions);

        assertEquals(605, results.size());
        assertEquals(2.0, first.ploidy(), EPSILON);
        assertEquals(1.0, first.normFactor(), EPSILON);
        assertEquals(1.0, first.diploidProportion(), EPSILON);
    }

    @Test
    public void twoRegions() throws Exception
    {
        List<FittingRegion> regions = Lists.newArrayList();
        regions.add(observedRegion(1, 1000, 100, 0.5, 1.0, 1.0));
        regions.add(observedRegion(1001, 2000, 100, 0.5, 1.0, 1.0));
        createFactoryAndRun(regions);

        assertEquals(605, results.size());
        assertEquals(2.0, first.ploidy(), EPSILON);
        assertEquals(1.0, first.normFactor(), EPSILON);
        assertEquals(1.0, first.diploidProportion(), EPSILON);
    }

    private FittingRegion observedRegion(int start, int end, int bafCount, double observedBAF, double observedTumorRatio, double observedNormalRatio)
    {
        return new ObservedRegion(
                chr1, start, end, true, SegmentSupport.NONE, bafCount, observedBAF, 10,
                observedTumorRatio, observedNormalRatio, observedNormalRatio, DIPLOID, false,
                0.5, 0, 0, 0, 0, 0,
                0, 2, 2, 0.5, 0, 0.5);
    }

    private void createFactoryAndRun(List<? extends FittingRegion> regions) throws Exception
    {
        factory = new FittedPurityFactory(mConfig, executorService, mCobaltChromosomes, mRegionFitCalculator, regions, List.of());
        factory.fitPurity();
        results = factory.getFittedPurities();
        var sorted = results.stream().sorted().toList();
        first = sorted.get(0);
        second = sorted.get(1);
        third = sorted.get(2);
    }
}
