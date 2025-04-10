package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.calcIndelInferredUnclippedPositions;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

import org.junit.Test;

public class ReadFiltersTest
{
    @Test
    public void testBasicFilters()
    {
        Junction posJunction = new Junction(CHR_1, 30, FORWARD);
        String readBases = REF_BASES_RANDOM_100.substring(10, 20); // irrelevant
        String readId = "READ_01";

        Read read = createRead(readId, 21, readBases, "10M10S");
        assertTrue(recordSoftClipsAndCrossesJunction(read, posJunction));

        read = createRead(readId, 23, readBases, "10M10S");
        assertTrue(recordSoftClipsAndCrossesJunction(read, posJunction));

        read = createRead(readId, 19, readBases, "10M10S");
        assertTrue(recordSoftClipsAndCrossesJunction(read, posJunction));

        // alignment too far from junction
        read = createRead(readId, 24, readBases, "10M10S");
        assertFalse(recordSoftClipsAndCrossesJunction(read, posJunction));

        // any realigned indel that crosses the junction is permitted
        read = createRead(readId, 10, readBases, "10M10I10M");
        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertTrue(recordSoftClipsAndCrossesJunction(read, posJunction));

        // a read with a low count of SNVs
        read = createRead(readId, 20, readBases, "12M");
        read.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(recordSoftClipsAndCrossesJunction(read, posJunction));

        Junction negJunction = new Junction(CHR_1, 30, REVERSE);

        read = createRead(readId, 30, readBases, "10S10M");
        assertTrue(recordSoftClipsAndCrossesJunction(read, negJunction));

        read = createRead(readId, 28, readBases, "10S10M");
        assertTrue(recordSoftClipsAndCrossesJunction(read, negJunction));

        read = createRead(readId, 32, readBases, "10S10M");
        assertTrue(recordSoftClipsAndCrossesJunction(read, negJunction));

        // alignment too far from junction
        read = createRead(readId, 33, readBases, "10S10M");
        assertFalse(recordSoftClipsAndCrossesJunction(read, negJunction));

        // any realigned indel that crosses the junction is permitted
        read = createRead(readId, 10, readBases, "10M10I10M");
        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertTrue(recordSoftClipsAndCrossesJunction(read, negJunction));

        // a read with a low count of SNVs
        read = createRead(readId, 28, readBases, "12M");
        read.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(recordSoftClipsAndCrossesJunction(read, negJunction));
    }
}
