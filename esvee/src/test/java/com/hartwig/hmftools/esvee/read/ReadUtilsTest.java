package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.read.Read.INVALID_INDEX;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ReadUtilsTest
{
    @Test
    public void testMiscReadFunctions()
    {
        String cigarStr = "5S20M8S";
        Read read = createSamRecord("READ_01", 20, REF_BASES.substring(15, 48), cigarStr);

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
    }
}
