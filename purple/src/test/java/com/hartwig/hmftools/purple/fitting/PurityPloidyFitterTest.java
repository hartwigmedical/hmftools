package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.FittingTestUtils.createObservedRegion;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Before;
import org.junit.Test;

public class PurityPloidyFitterTest extends FittingTestBase
{
    @Before
    public void setup()
    {
        super.setup();
    }

    private List<ObservedRegion> buildDefaultObservedRegions()
    {
        List<ObservedRegion> observedRegions = Lists.newArrayList();

        Chromosome chr1 = HumanChromosome._1;

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1, 1000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1001, 2000, 0.5, 0.5, GermlineStatus.DIPLOID, 2));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 2001, 3000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 3001, 4000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        return observedRegions;
    }

    @Test
    public void testSomaticFit()
    {
        List<ObservedRegion> observedRegions = buildDefaultObservedRegions();

        PurityPloidyFitter fitter = new PurityPloidyFitter(
                mConfig, mReferenceData, mSampleData, null, mRegionFitCalculator, observedRegions, Gender.FEMALE, false);

        assertTrue(fitter.isValid());

        fitter.run();

        assertNotNull(fitter.finalFit());
    }
}
