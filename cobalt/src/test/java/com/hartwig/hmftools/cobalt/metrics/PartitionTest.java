package com.hartwig.hmftools.cobalt.metrics;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.cobalt.ChromosomeData;

import org.junit.Test;

public class PartitionTest
{
    @Test
    public void createPartitionsTest()
    {
        Partition partition = new Partition("21", 2_000_001, 3_000_000, 10_000);

        assertEquals(100, partition.TargetRegions.size());

        TargetRegionData region0 = partition.TargetRegions.get(0);
        assertEquals(2_000_001, region0.start());
        assertEquals(2_010_000, region0.end());

        TargetRegionData region1 = partition.TargetRegions.get(1);
        assertEquals(2_010_001, region1.start());
        assertEquals(2_020_000, region1.end());
    }
}
