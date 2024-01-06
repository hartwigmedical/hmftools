package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.old.PrimaryAssembler.realignForJunction;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class PrimaryAssemblyTest
{
    @Test
    public void testReadIndelToSoftClips()
    {
        Junction fwdJunction = new Junction(CHR_1, 30, POS_STRAND);

        Read read = createSamRecord(TEST_READ_ID, 20, REF_BASES.substring(20, 40), "10M2I8M");
        assertEquals(12, read.getReadIndexAtReferencePosition(fwdJunction.Position)); // being the first base after the insert

        Read realignedRead = realignForJunction(read, fwdJunction);
        assertEquals(10, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(29, realignedRead.getAlignmentEnd());
        assertEquals(39, realignedRead.getUnclippedEnd());

        Junction revJunction = new Junction(CHR_1, 30, NEG_STRAND);

        read = createSamRecord(
                TEST_READ_ID, 23, REF_BASES.substring(23, 31) + "AA" + REF_BASES.substring(31, 40), "8M2I10M");

        assertEquals(7, read.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        realignedRead = realignForJunction(read, revJunction);
        assertEquals(9, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(31, realignedRead.getAlignmentStart());
        assertEquals(31, realignedRead.getAlignmentStart());
        assertEquals(21, realignedRead.getUnclippedStart());
    }



}
