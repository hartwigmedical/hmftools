package com.hartwig.hmftools.cobalt.lowcov;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.junit.Test;

public class LowCoverageRatioBuilderTest
{
    @Test
    public void testCalcConsoliationCount()
    {
        assertEquals(1, LowCoverageRatioBuilder.calcConsolidationCount(100.0));
        assertEquals(1, LowCoverageRatioBuilder.calcConsolidationCount(51.0));

        // starting from 50 reads we consolidate
        assertEquals(10, LowCoverageRatioBuilder.calcConsolidationCount(50.0));
        assertEquals(30, LowCoverageRatioBuilder.calcConsolidationCount(20.0));
        assertEquals(100, LowCoverageRatioBuilder.calcConsolidationCount(5.0));
        assertEquals(300, LowCoverageRatioBuilder.calcConsolidationCount(2.0));
        assertEquals(500, LowCoverageRatioBuilder.calcConsolidationCount(1.0));
        assertEquals(1000, LowCoverageRatioBuilder.calcConsolidationCount(0.5));
        assertEquals(1000, LowCoverageRatioBuilder.calcConsolidationCount(0.05));
        assertEquals(1000, LowCoverageRatioBuilder.calcConsolidationCount(0.001));
    }

    @Test
    public void testCalcConsolidateBoundaries()
    {
        List<Integer> windowPositions = new ArrayList<>();
        windowPositions.add(1_001);
        windowPositions.add(2_001);
        windowPositions.add(5_001);
        windowPositions.add(6_001);
        windowPositions.add(10_001);

        // test we do not span over centromere
        windowPositions.add(3_020_001);
        windowPositions.add(3_021_001);
        windowPositions.add(3_024_001);
        windowPositions.add(3_029_001);

        List<LowCovBucket> buckets = LowCoverageRatioBuilder.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).startPosition);
        assertEquals(8_001, buckets.get(0).endPosition);
        assertEquals(9_001, buckets.get(1).startPosition);
        assertEquals(11_001, buckets.get(1).endPosition);
        assertEquals(3_020_001, buckets.get(2).startPosition);
        assertEquals(3_030_001, buckets.get(2).endPosition);

        windowPositions.add(3_034_001);

        buckets = LowCoverageRatioBuilder.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(4, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).startPosition);
        assertEquals(8_001, buckets.get(0).endPosition);
        assertEquals(9_001, buckets.get(1).startPosition);
        assertEquals(11_001, buckets.get(1).endPosition);
        assertEquals(3_020_001, buckets.get(2).startPosition);
        assertEquals(3_031_001, buckets.get(2).endPosition);
        assertEquals(3_032_001, buckets.get(3).startPosition);
        assertEquals(3_035_001, buckets.get(3).endPosition);

    }

    @Test
    public void testCalcConsolidateBoundaryRatios()
    {
        final ListMultimap<Chromosome, ReadRatio> rawRatios = ArrayListMultimap.create();

        Chromosome chr = new Chromosome("chr1", 23310);


        // add in some chromosome read ratio
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(1001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(2001).ratio(-1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(3001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(5001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(9001).ratio(1.0).build());

        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(10001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(12001).ratio(-1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(13001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(14001).ratio(1.0).build());
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(16001).ratio(1.0).build());

        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(19001).ratio(1.0).build());

        List<LowCovBucket> buckets = Objects.requireNonNull(LowCoverageRatioBuilder.consolidateIntoBuckets(rawRatios, 4)).get(chr);

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1001, buckets.get(0).startPosition);
        assertEquals(5001, buckets.get(0).bucketPosition);
        assertEquals(9001, buckets.get(0).endPosition);
        assertEquals(10001, buckets.get(1).startPosition);
        assertEquals(13001, buckets.get(1).bucketPosition);
        assertEquals(17001, buckets.get(1).endPosition);
        assertEquals(18001, buckets.get(2).startPosition);
        assertEquals(19001, buckets.get(2).bucketPosition);
        assertEquals(20001, buckets.get(2).endPosition);

        // put a masked out ratio at the end, should also work
        rawRatios.put(chr, ImmutableReadRatio.builder().chromosome(chr.contig).position(20001).ratio(-1.0).build());

        buckets = Objects.requireNonNull(LowCoverageRatioBuilder.consolidateIntoBuckets(rawRatios, 4)).get(chr);

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1001, buckets.get(0).startPosition);
        assertEquals(9001, buckets.get(0).endPosition);
        assertEquals(10001, buckets.get(1).startPosition);
        assertEquals(17001, buckets.get(1).endPosition);
        assertEquals(18001, buckets.get(2).startPosition);
        assertEquals(20001, buckets.get(2).endPosition);
    }

}
