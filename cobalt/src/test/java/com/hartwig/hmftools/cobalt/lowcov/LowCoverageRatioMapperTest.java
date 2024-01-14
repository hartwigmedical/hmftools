package com.hartwig.hmftools.cobalt.lowcov;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.cobalt.CobaltColumns;

import org.junit.Test;

import tech.tablesaw.api.*;

public class LowCoverageRatioMapperTest
{
    @Test
    public void testCalcConsoliationCount()
    {
        assertEquals(1, LowCoverageRatioMapper.calcConsolidationCount(16.0));
        assertEquals(1, LowCoverageRatioMapper.calcConsolidationCount(8.1));

        // starting from 8.0 reads depth we consolidate
        assertEquals(10, LowCoverageRatioMapper.calcConsolidationCount(8.0));
        assertEquals(30, LowCoverageRatioMapper.calcConsolidationCount(3.0));
        assertEquals(100, LowCoverageRatioMapper.calcConsolidationCount(0.8));
        assertEquals(300, LowCoverageRatioMapper.calcConsolidationCount(0.3));
        assertEquals(500, LowCoverageRatioMapper.calcConsolidationCount(0.15));
        assertEquals(1000, LowCoverageRatioMapper.calcConsolidationCount(0.08));
        assertEquals(1000, LowCoverageRatioMapper.calcConsolidationCount(0.008));
        assertEquals(1000, LowCoverageRatioMapper.calcConsolidationCount(0.00015));
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

        List<LowCovBucket> buckets = LowCoverageRatioMapper.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).StartPosition);
        assertEquals(8_001, buckets.get(0).EndPosition);
        assertEquals(9_001, buckets.get(1).StartPosition);
        assertEquals(11_001, buckets.get(1).EndPosition);
        assertEquals(3_020_001, buckets.get(2).StartPosition);
        assertEquals(3_030_001, buckets.get(2).EndPosition);

        windowPositions.add(3_034_001);

        buckets = LowCoverageRatioMapper.consolidateIntoBuckets(windowPositions, 4);

        assertEquals(4, buckets.size());

        // check each one
        assertEquals(1_001, buckets.get(0).StartPosition);
        assertEquals(8_001, buckets.get(0).EndPosition);
        assertEquals(9_001, buckets.get(1).StartPosition);
        assertEquals(11_001, buckets.get(1).EndPosition);
        assertEquals(3_020_001, buckets.get(2).StartPosition);
        assertEquals(3_031_001, buckets.get(2).EndPosition);
        assertEquals(3_032_001, buckets.get(3).StartPosition);
        assertEquals(3_035_001, buckets.get(3).EndPosition);

    }

    @Test
    public void testCalcConsolidateBoundaryRatios()
    {
        final Table rawRatios = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.RATIO));

        Row row;
        // add in some chromosome read ratio
        appendReadRatio(rawRatios, "chr1", 1001, 1.0);
        appendReadRatio(rawRatios, "chr1", 2001, -1.0);
        appendReadRatio(rawRatios, "chr1", 3001, 1.0);
        appendReadRatio(rawRatios, "chr1", 5001, 1.0);
        appendReadRatio(rawRatios, "chr1", 9001, 1.0);

        appendReadRatio(rawRatios, "chr1", 10001, 1.0);
        appendReadRatio(rawRatios, "chr1", 12001, -1.0);
        appendReadRatio(rawRatios, "chr1", 13001, 1.0);
        appendReadRatio(rawRatios, "chr1", 14001, 1.0);
        appendReadRatio(rawRatios, "chr1", 16001, 1.0);

        appendReadRatio(rawRatios, "chr1", 19001, 1.0);

        List<LowCovBucket> buckets = Objects.requireNonNull(LowCoverageRatioMapper.consolidateIntoBuckets(rawRatios, 4)).get("chr1");

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1001, buckets.get(0).StartPosition);
        assertEquals(5001, buckets.get(0).BucketPosition);
        assertEquals(9001, buckets.get(0).EndPosition);
        assertEquals(10001, buckets.get(1).StartPosition);
        assertEquals(13001, buckets.get(1).BucketPosition);
        assertEquals(17001, buckets.get(1).EndPosition);
        assertEquals(18001, buckets.get(2).StartPosition);
        assertEquals(19001, buckets.get(2).BucketPosition);
        assertEquals(20001, buckets.get(2).EndPosition);

        // put a masked out ratio at the end, should also work
        row = rawRatios.appendRow(); row.setString("chromosome", "chr1"); row.setInt("position", 20001); row.setDouble("ratio", -1.0);

        buckets = Objects.requireNonNull(LowCoverageRatioMapper.consolidateIntoBuckets(rawRatios, 4)).get("chr1");

        assertEquals(3, buckets.size());

        // check each one
        assertEquals(1001, buckets.get(0).StartPosition);
        assertEquals(9001, buckets.get(0).EndPosition);
        assertEquals(10001, buckets.get(1).StartPosition);
        assertEquals(17001, buckets.get(1).EndPosition);
        assertEquals(18001, buckets.get(2).StartPosition);
        assertEquals(20001, buckets.get(2).EndPosition);
    }

    @SuppressWarnings("SameParameterValue")
    private static void appendReadRatio(Table table, String chromosome, int position, double ratio)
    {
        Row row = table.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, chromosome);
        row.setInt(CobaltColumns.POSITION, position);
        row.setDouble(CobaltColumns.RATIO, ratio);
    }
}
