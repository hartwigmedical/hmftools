package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.read.ReadFilters.recordSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.types.Junction;

import org.junit.Test;

public class ReadFiltersTest
{
    @Test
    public void testBasicFilters()
    {
        Read read = createRead("READ_01", 20, REF_BASES_RANDOM_100.substring(15, 48), "5S20M1S");

        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 19, POS_STRAND)));
        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 22, NEG_STRAND)));
        assertFalse(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 15, NEG_STRAND))); // doesn't extend far enough

        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 39, NEG_STRAND)));
        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 39, POS_STRAND)));
        assertFalse(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 40, POS_STRAND))); // doesn't extend far enough
    }
}
