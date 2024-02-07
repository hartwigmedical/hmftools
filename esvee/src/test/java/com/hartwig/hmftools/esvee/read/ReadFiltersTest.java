package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.read.ReadFilters.recordSoftClipsNearJunction;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.common.Junction;

import org.junit.Test;

public class ReadFiltersTest
{
    @Test
    public void testBasicFilters()
    {
        Read read = createSamRecord("READ_01", 20, REF_BASES.substring(15, 48), "5S20M1S");

        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 19, POS_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 18, NEG_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 22, NEG_STRAND)));
        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 17, NEG_STRAND)));

        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 39, NEG_STRAND)));
        assertTrue(recordSoftClipsNearJunction(read, new Junction(CHR_1, 39, POS_STRAND)));
        assertFalse(recordSoftClipsNearJunction(read, new Junction(CHR_1, 40, POS_STRAND))); // doesn't extend far enough
    }
}
