package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class ReadUtilsTest
{
    @Test
    public void testMiscReadFunctions()
    {
        String cigarStr = "5S20M8S";
        Read read = createRead(TEST_READ_ID, 20, REF_BASES_RANDOM_100.substring(15, 48), cigarStr);

        assertEquals(15, read.unclippedStart());
        assertEquals(47, read.unclippedEnd());
        assertEquals(3, read.cigarElements().size());
        assertEquals(cigarStr, read.cigarString());

        assertEquals(5, read.getReadIndexAtReferencePosition(20));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(19));
        assertEquals(4, read.getReadIndexAtReferencePosition(19, true));
        assertEquals(0, read.getReadIndexAtReferencePosition(15, true));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(14, true));

        assertEquals(24, read.getReadIndexAtReferencePosition(39));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(40));
        assertEquals(25, read.getReadIndexAtReferencePosition(40, true));
        assertEquals(32, read.getReadIndexAtReferencePosition(47, true));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(48, true));

        // with INDELs

        cigarStr = "4S20M10D10M5I10M3S";
        // positions 4S then 20-39 = M,   40-49 = D, 50-59 = M,    T at 59,     60-69 = M, then 3S
        // index: 4S is 0-3, 20M is 4-23, D is none, 10M is 24-33, 5I is 34-38, 10M is 39-48, 3S is 49-51
        read = createRead(
                TEST_READ_ID, 20,
                REF_BASES_RANDOM_100.substring(16, 40) + REF_BASES_RANDOM_100.substring(50, 60) + "GGGGG" + REF_BASES_RANDOM_100.substring(60, 73),
                cigarStr);

        assertEquals(16, read.unclippedStart());
        assertEquals(69, read.alignmentEnd());
        assertEquals(72, read.unclippedEnd());
        assertEquals(7, read.cigarElements().size());
        assertEquals(52, read.basesLength());

        assertEquals(0, read.getReadIndexAtReferencePosition(16, true));
        assertEquals(4, read.getReadIndexAtReferencePosition(20, false));
        assertEquals(5, read.getReadIndexAtReferencePosition(21, false));
        assertEquals(23, read.getReadIndexAtReferencePosition(39, false)); // start of DEL
        assertEquals(23, read.getReadIndexAtReferencePosition(48, false)); // within DEL

        assertEquals(24, read.getReadIndexAtReferencePosition(50, false)); // end of DEL
        assertEquals(33, read.getReadIndexAtReferencePosition(59, false)); // end of next M
        assertEquals(39, read.getReadIndexAtReferencePosition(60, false)); // base after insert
        assertEquals(48, read.getReadIndexAtReferencePosition(69, false)); // final aligned base
        assertEquals(49, read.getReadIndexAtReferencePosition(70, true));
        assertEquals(50, read.getReadIndexAtReferencePosition(71, true));
        assertEquals(51, read.getReadIndexAtReferencePosition(72, true)); // final SC base

        // compare to SAM record method, which is 1-based
        assertEquals(5, read.bamRecord().getReadPositionAtReferencePosition(20));
        assertEquals(24, read.bamRecord().getReadPositionAtReferencePosition(39)); // start of DEL

        assertEquals(25, read.bamRecord().getReadPositionAtReferencePosition(50)); // end of DEL
        assertEquals(34, read.bamRecord().getReadPositionAtReferencePosition(59)); // end of next M
        assertEquals(40, read.bamRecord().getReadPositionAtReferencePosition(60)); // base after insert
        assertEquals(49, read.bamRecord().getReadPositionAtReferencePosition(69)); // final aligned base
    }
}
