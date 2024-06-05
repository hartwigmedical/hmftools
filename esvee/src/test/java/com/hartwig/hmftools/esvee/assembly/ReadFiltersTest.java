package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

import org.junit.Test;

public class ReadFiltersTest
{
    @Test
    public void testBasicFilters()
    {
        Read read = createRead("READ_01", 20, REF_BASES_RANDOM_100.substring(15, 48), "5S20M1S");

        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 19, FORWARD)));
        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 22, REVERSE)));
        assertFalse(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 15, REVERSE))); // doesn't extend far enough

        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 39, REVERSE)));
        assertTrue(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 39, FORWARD)));
        assertFalse(recordSoftClipsAndCrossesJunction(read, new Junction(CHR_1, 40, FORWARD))); // doesn't extend far enough
    }
}
