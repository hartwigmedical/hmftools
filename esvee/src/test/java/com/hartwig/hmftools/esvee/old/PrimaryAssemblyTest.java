package com.hartwig.hmftools.esvee.old;

import static org.junit.Assert.assertEquals;

public class PrimaryAssemblyTest
{
    /*
    @Test
    public void testReadIndelToSoftClips()
    {
        Junction fwdJunction = new Junction(CHR_1, 30, POS_STRAND);

        Read read = createSamRecord(TEST_READ_ID, 20, REF_BASES.substring(20, 40), "10M2I8M");
        assertEquals(12, read.getReadIndexAtReferencePosition(fwdJunction.Position)); // being the first base after the insert

        Read realignedRead = realignForJunction(read, fwdJunction);
        assertEquals(10, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(29, realignedRead.alignmentEnd());
        assertEquals(39, realignedRead.unclippedEnd());

        Junction revJunction = new Junction(CHR_1, 30, NEG_STRAND);

        read = createSamRecord(
                TEST_READ_ID, 23, REF_BASES.substring(23, 31) + "AA" + REF_BASES.substring(31, 40), "8M2I10M");

        assertEquals(7, read.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        realignedRead = realignForJunction(read, revJunction);
        assertEquals(9, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(31, realignedRead.alignmentStart());
        assertEquals(31, realignedRead.alignmentStart());
        assertEquals(21, realignedRead.unclippedStart());
    }
    */
}
