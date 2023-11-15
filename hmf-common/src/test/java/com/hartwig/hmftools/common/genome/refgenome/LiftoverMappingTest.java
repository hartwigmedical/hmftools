package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

public class LiftoverMappingTest
{
    @Test
    public void testConversions()
    {
        final InputStream inputStream = LiftoverMappingTest.class.getResourceAsStream("/genome/liftover_mappings.tsv");
        List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());

        GenomeLiftoverCache genomeLiftoverCache = new GenomeLiftoverCache();
        genomeLiftoverCache.loadFile(lines);
        assertTrue(genomeLiftoverCache.hasMappings());

        assertEquals(UNMAPPED_POSITION, genomeLiftoverCache.convertPosition(CHR_1, 5000));

        assertEquals(10001, genomeLiftoverCache.convertPosition(CHR_1, 10001));

        // test a reversed region
        assertEquals(501617, genomeLiftoverCache.convertPosition(CHR_1, 317720));
        assertEquals(347969, genomeLiftoverCache.convertPosition(CHR_1, 471368));

        // test reversing these positions
        assertEquals(317720, genomeLiftoverCache.convertPositionTo37(CHR_1, 501617));
        assertEquals(471368, genomeLiftoverCache.convertPositionTo37(CHR_1, 347969));
        assertEquals(10001, genomeLiftoverCache.convertPositionTo37(CHR_1, 10001));
    }
}
