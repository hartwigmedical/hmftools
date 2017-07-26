package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.pcf.ImmutablePCFRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.junit.Before;
import org.junit.Test;

public class SegmentMergeTest {

    private List<GenomeRegion> normal;
    private List<GenomeRegion> tumor;

    @Before
    public void setup() {
        normal = Lists.newArrayList();
        tumor = Lists.newArrayList();
    }

    @Test
    public void testNoNewRegions() {
        GenomeRegion normalRegion1 = create("1", 1, 1000);
        GenomeRegion normalRegion2 = create("1", 1001, 2000);
        GenomeRegion normalRegion3 = create("1", 2001, 3000);
        GenomeRegion tumorRegion1 = create("1", 1001, 2000);
        GenomeRegion tumorRegion2 = create("1", 2001, 3000);
        GenomeRegion tumorRegion3 = create("2", 1001, 2000);

        normal.add(normalRegion1);
        normal.add(normalRegion2);
        normal.add(normalRegion3);
        tumor.add(tumorRegion1);
        tumor.add(tumorRegion2);
        tumor.add(tumorRegion3);
        List<GenomeRegion> regions = SegmentMerge.merge(normal, tumor);
        assertEquals(4, regions.size());
        assertEquals(normalRegion1, regions.get(0));
        assertEquals(normalRegion2, regions.get(1));
        assertEquals(normalRegion3, regions.get(2));
        assertEquals(tumorRegion3, regions.get(3));
    }

    @Test
    public void testOverlaps() {
        GenomeRegion normalRegion1 = create("1", 1, 10000);
        GenomeRegion tumorRegion1 = create("1", 1001, 2000);
        GenomeRegion tumorRegion2 = create("1", 2001, 3000);
        GenomeRegion tumorRegion3 = create("1", 5001, 6000);
        GenomeRegion tumorRegion4 = create("1", 9001, 11000);

        normal.add(normalRegion1);
        tumor.add(tumorRegion1);
        tumor.add(tumorRegion2);
        tumor.add(tumorRegion3);
        tumor.add(tumorRegion4);
        List<GenomeRegion> regions = SegmentMerge.merge(normal, tumor);
        assertEquals(8, regions.size());
        assertEquals(create("1", 1, 1000), regions.get(0));
        assertEquals(tumorRegion1, regions.get(1));
        assertEquals(tumorRegion2, regions.get(2));
        assertEquals(create("1", 3001, 5000), regions.get(3));
        assertEquals(tumorRegion3, regions.get(4));
        assertEquals(create("1", 6001, 9000), regions.get(5));
        assertEquals(create("1", 9001, 10000), regions.get(6));
        assertEquals(create("1", 10001, 11000), regions.get(7));
    }

    private GenomeRegion create(String chromosome, long start, long end) {
        return ImmutablePCFRegion.builder().chromosome(chromosome).start(start).end(end).build();
    }
}
