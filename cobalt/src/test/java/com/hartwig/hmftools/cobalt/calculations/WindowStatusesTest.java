package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class WindowStatusesTest extends CalculationsTestBase
{

    @Test
    public void excludeTest()
    {
        ListMultimap<Chromosome, GCProfile> gcProfiles = ArrayListMultimap.create();
        gcProfiles.put(_1, gcProfile(_1, 1, 0.47, 0.50));
        gcProfiles.put(_1, gcProfile(_1, 1001, 0.47, 0.50));
        gcProfiles.put(_1, gcProfile(_1, 2001, 0.47, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 3001, 0.47, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 4001, 0.47, 1.0));
        gcProfiles.put(_1, gcProfile(_1, 5001, 0.47, 1.0));

        gcProfiles.put(_2, gcProfile(_2, 1, 0.57, 1.0));
        gcProfiles.put(_2, gcProfile(_2, 1001, 0.57, 1.0));
        gcProfiles.put(_2, gcProfile(_2, 2001, 0.57, 0.50));
        gcProfiles.put(_2, gcProfile(_2, 3001, 0.57, 1.0));

        List<ChrBaseRegion> exclusions = new ArrayList<>();
        exclusions.add(cbr(_1, 1100, 1200));
        exclusions.add(cbr(_1, 1300, 1400));
        exclusions.add(cbr(_1, 1100, 1200));
        exclusions.add(cbr(_1, 4100, 4200));

        WindowStatuses statuses = new WindowStatuses(gcProfiles,exclusions);
        assertTrue(statuses.exclude(_1, dr(_1, 1001, 0.47, 0.50)));
        assertFalse(statuses.exclude(_1, dr(_1, 2001, 0.47, 0.50)));
        assertFalse(statuses.exclude(_1, dr(_1, 3001, 0.47, 0.50)));
        assertTrue(statuses.exclude(_1, dr(_1, 4001, 0.47, 0.50)));
        assertFalse(statuses.exclude(_1, dr(_1, 5001, 0.47, 0.50)));
        assertFalse(statuses.exclude(_2, dr(_2, 1001, 0.47, 0.50)));
        assertTrue(statuses.exclude(_2, dr(_2, 2001, 0.47, 0.50)));
        assertFalse(statuses.exclude(_2, dr(_2, 3001, 0.47, 0.50)));
    }
}
