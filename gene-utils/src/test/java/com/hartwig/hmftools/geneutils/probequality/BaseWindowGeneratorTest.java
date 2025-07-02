package com.hartwig.hmftools.geneutils.probequality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class BaseWindowGeneratorTest
{
    private final BaseWindowGenerator mGen;

    public BaseWindowGeneratorTest()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.putAll(Map.of(
                "10", "AACCGGTTAACCGGTT",
                "15", "TTGGCCAATTGGCCAATTGGCCAAT",
                "17", "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
        ));
        // Compute chromosome lengths based on base sequences.
        refGenome.ChromosomeLengths.putAll(
                refGenome.RefGenomeMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().length())));

        SpecificRegions specificRegions = new SpecificRegions();
        // Include all of chr 10
        specificRegions.Chromosomes.add("10");
        // Include part of chr 15
        specificRegions.addRegion(new ChrBaseRegion("15", 3, 8));
        specificRegions.addRegion(new ChrBaseRegion("15", 5, 9));   // Check overlap behaviour
        specificRegions.addRegion(new ChrBaseRegion("15", 20, 24));
        // Exclude all of chr 17

        // Windows will be: [1, 5], [4, 8], [7, 11], [10, 14], [13, 17], [16, 20], [19, 23], [22, 26], ...
        int baseWindowLength = 5;
        int baseWindowSpacing = 3;
        int batchSize = 3;

        mGen = new BaseWindowGenerator(refGenome, specificRegions, baseWindowLength, baseWindowSpacing, batchSize);
    }

    @Test
    public void test_createBaseWindowBatches() {
        List<List<BaseWindowGenerator.BaseWindow>> expected = List.of(
                List.of(
                    new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("10", 1, 5), "AACCG".getBytes()),
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("10", 4, 8), "CGGTT".getBytes()),
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("10", 7, 11), "TTAAC".getBytes())
                ),
                List.of(
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("10", 10, 14), "ACCGG".getBytes()),
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("15", 1, 5), "TTGGC".getBytes()),
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("15", 4, 8), "GCCAA".getBytes())
                        ),
                List.of(
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("15", 7, 11), "AATTG".getBytes()),
                        new BaseWindowGenerator.BaseWindow(new ChrBaseRegion("15", 19, 23), "GGCCA".getBytes())
                )
        );
        List<List<BaseWindowGenerator.BaseWindow>> actual = mGen.createBaseWindowBatches().toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_chromosome_all()
    {
        String chr = "10";
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion(chr, 1, 5),
                new ChrBaseRegion(chr, 4, 8),
                new ChrBaseRegion(chr, 7, 11),
                new ChrBaseRegion(chr, 10, 14)
                // Note last bit of chromosome is excluded since the window would extend past the end.
        );
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions(chr).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_chromosome_part()
    {
        String chr = "15";
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion(chr, 1, 5),
                new ChrBaseRegion(chr, 4, 8),
                new ChrBaseRegion(chr, 7, 11),
                new ChrBaseRegion(chr, 19, 23)
                // Last window excluded due to extending past chromosome end.
        );
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions(chr).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_chromosome_excluded()
    {
        List<ChrBaseRegion> expected = List.of();
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions("17").toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_region_empty()
    {
        ChrBaseRegion region = new ChrBaseRegion("10", 11, 10);
        List<ChrBaseRegion> expected = List.of();
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions(region).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_region_single()
    {
        String chr = "10";
        ChrBaseRegion region = new ChrBaseRegion(chr, 11, 12);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion(chr, 10, 14)
        );
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions(region).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_createBaseWindowRegions_region_multiple()
    {
        String chr = "10";
        ChrBaseRegion region = new ChrBaseRegion(chr, 12, 30);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion(chr, 10, 14),
                new ChrBaseRegion(chr, 13, 17),
                new ChrBaseRegion(chr, 16, 20),
                new ChrBaseRegion(chr, 19, 23),
                new ChrBaseRegion(chr, 22, 26),
                new ChrBaseRegion(chr, 25, 29),
                new ChrBaseRegion(chr, 28, 32)
        );
        List<ChrBaseRegion> actual = mGen.createBaseWindowRegions(region).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_baseWindowStartCoveringPosition()
    {
        assertEquals(1, mGen.baseWindowStartCoveringPosition(1));
        assertEquals(1, mGen.baseWindowStartCoveringPosition(2));
        assertEquals(1, mGen.baseWindowStartCoveringPosition(3));
        assertEquals(4, mGen.baseWindowStartCoveringPosition(4));
        assertEquals(4, mGen.baseWindowStartCoveringPosition(5));
        assertEquals(4, mGen.baseWindowStartCoveringPosition(6));
        assertEquals(7, mGen.baseWindowStartCoveringPosition(8));
        assertEquals(7, mGen.baseWindowStartCoveringPosition(9));
        assertEquals(10, mGen.baseWindowStartCoveringPosition(10));
        assertEquals(13, mGen.baseWindowStartCoveringPosition(14));
    }

    @Test
    public void test_batchRegions_single()
    {
        List<ChrBaseRegion> regions = List.of(
                new ChrBaseRegion("1", 1, 10),
                new ChrBaseRegion("2", 10, 20)
        );
        List<List<ChrBaseRegion>> expected = List.of(
                List.of(regions.get(0), regions.get(1))
        );
        List<List<ChrBaseRegion>> actual = mGen.batchRegions(regions.stream()).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_batchRegions_multiple()
    {
        List<ChrBaseRegion> regions = List.of(
                new ChrBaseRegion("1", 1, 10),
                new ChrBaseRegion("2", 10, 20),
                new ChrBaseRegion("3", 20, 30),
                new ChrBaseRegion("4", 30, 40),
                new ChrBaseRegion("5", 40, 50),
                new ChrBaseRegion("6", 50, 60),
                new ChrBaseRegion("7", 60, 70)
        );
        List<List<ChrBaseRegion>> expected = List.of(
                List.of(regions.get(0), regions.get(1), regions.get(2)),
                List.of(regions.get(3), regions.get(4), regions.get(5)),
                List.of(regions.get(6))
        );
        List<List<ChrBaseRegion>> actual = mGen.batchRegions(regions.stream()).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_isSequenceNormal()
    {
        assertTrue(BaseWindowGenerator.isSequenceNormal("AAAAGCTaGtCAGTgcgtacg".getBytes()));
        assertFalse(BaseWindowGenerator.isSequenceNormal("GGGGGNNNNNNNNNNTCTC".getBytes()));
        assertFalse(BaseWindowGenerator.isSequenceNormal("GCTGUTCAGCTxxxxTAC".getBytes()));
    }
}
