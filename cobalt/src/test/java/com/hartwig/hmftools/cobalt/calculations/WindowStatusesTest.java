package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.diploid.DiploidStatus;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Before;
import org.junit.Test;

public class WindowStatusesTest extends CalculationsTestBase
{
    private WindowStatuses statuses;

    @Before
    public void setup()
    {
        ListMultimap<Chromosome, GCProfile> gcProfiles = ArrayListMultimap.create();
        gcProfiles.put(_1, gcProfile(_1, 1, 0.47, 0.70));    // not mappable
        gcProfiles.put(_1, gcProfile(_1, 1001, 0.48, 0.50)); // not mappable
        gcProfiles.put(_1, gcProfile(_1, 2001, 0.49, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 3001, 0.50, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 4001, 0.51, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 5001, 0.52, 1.0));

        gcProfiles.put(_2, gcProfile(_2, 1, 0.57, 1.0));
        gcProfiles.put(_2, gcProfile(_2, 1001, 0.58, 1.0));
        gcProfiles.put(_2, gcProfile(_2, 2001, 0.59, 0.50)); // not mappable
        gcProfiles.put(_2, gcProfile(_2, 3001, 0.57, 1.0));

        gcProfiles.put(_3, gcProfile(_3, 1, 0.47, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 1001, 0.48, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 2001, 0.49, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 3001, 0.50, 0.8)); // not mappable
        gcProfiles.put(_3, gcProfile(_3, 4001, 0.47, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 5001, 0.48, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 6001, 0.49, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 7001, 0.50, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 8001, 0.47, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 9001, 0.48, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 10001, 0.49, 1.0));
        gcProfiles.put(_3, gcProfile(_3, 11001, 0.50, 1.0));

        List<ChrBaseRegion> exclusions = new ArrayList<>();
        exclusions.add(cbr(_1, 1100, 1200));
        exclusions.add(cbr(_1, 1300, 1400));
        exclusions.add(cbr(_1, 4100, 4200));
        exclusions.add(cbr(_3, 1100, 1200));
        exclusions.add(cbr(_3, 7100, 7200));

        ListMultimap<Chromosome, DiploidStatus> diploidRegions = ArrayListMultimap.create();
        for(int i = 0; i < 7; i++)
        {
            diploidRegions.put(_1, ds(_1, 1000 * i + 1, true));
        }
        for(int i = 0; i < 4; i++)
        {
            diploidRegions.put(_2, ds(_2, 1000 * i + 1, true));
        }
        diploidRegions.put(_3, ds(_3, 1, true));
        diploidRegions.put(_3, ds(_3, 1001, true));
        diploidRegions.put(_3, ds(_3, 2001, true));
        diploidRegions.put(_3, ds(_3, 3001, false));
        diploidRegions.put(_3, ds(_3, 4001, true));
        diploidRegions.put(_3, ds(_3, 5001, true));
        diploidRegions.put(_3, ds(_3, 6001, false));
        diploidRegions.put(_3, ds(_3, 7001, false));
        diploidRegions.put(_3, ds(_3, 8001, true));
        diploidRegions.put(_3, ds(_3, 9001, true));
        diploidRegions.put(_3, ds(_3, 10001, true));

        statuses = new WindowStatuses(gcProfiles, exclusions, diploidRegions);
    }

    @Test
    public void gcReferenceValueTest()
    {
        assertEquals(0.47, statuses.referenceGcValueForWindow(_1, 0), 0.0001);
        assertEquals(0.48, statuses.referenceGcValueForWindow(_1, 1000), 0.0001);
        assertEquals(0.57, statuses.referenceGcValueForWindow(_2, 1), 0.0001);
        assertEquals(0.47, statuses.referenceGcValueForWindow(_3, 1), 0.0001);
        assertEquals(0.48, statuses.referenceGcValueForWindow(_3, 1000), 0.0001);
    }

    @Test
    public void excludeTest()
    {
        checkExcluded(_1, 1);    // not mappable, contains excluded regions
        checkExcluded(_1, 1001); // not mappable
        checkIncluded(_1, 2001);
        checkIncluded(_1, 3001);
        checkExcluded(_1, 4001); // contains excluded regions
        checkIncluded(_1, 5001);

        checkIncluded(_2, 1);
        checkIncluded(_2, 1001);
        checkExcluded(_2, 2001); // not mappable
        checkIncluded(_2, 3001);

        checkIncluded(_3, 1);
        checkExcluded(_3, 1001); // contains excluded regions
        checkIncluded(_3, 2001);
        checkExcluded(_3, 3001); // not mappable, not diploid
        checkIncluded(_3, 4001);
        checkIncluded(_3, 5001);
        checkExcluded(_3, 6001); // not diploid
        checkExcluded(_3, 7001); // contains excluded regions, not diploid
        checkIncluded(_3, 8001);
        checkIncluded(_3, 9001);
        checkIncluded(_3, 10001);
        checkExcluded(_3, 11001); // beyond the range of the diploid regions data
    }

    private DiploidStatus ds(Chromosome chromosome, int position, boolean status)
    {
        return new DiploidStatus(chromosome.contig(), position, position + 999, status);
    }

    private void checkExcluded(Chromosome chromosome, int position)
    {
        assertTrue(statuses.exclude(chromosome, dr(chromosome, position)));
    }

    private void checkIncluded(Chromosome chromosome, int position)
    {
        assertFalse(statuses.exclude(chromosome, dr(chromosome, position)));
    }

    private DepthReading dr(Chromosome chromosome, int position)
    {
        return dr(chromosome, position, 0.47, 0.50);
    }
}
