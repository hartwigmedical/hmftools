package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndelIndices;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class SbxBamUtilsTest
{
    @Test
    public void
    testGetDuplexIndelsNoDuplexIndels()
    {
        String ycTagStr = "+100+"; // no leading simplex count or trailing tail, both meaning zero
        List<Integer> duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertTrue(duplexIndelIndices.isEmpty());

        ycTagStr = "10+100I100+"; // no trailing tail, meaning zero
        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(110));

        ycTagStr = "+100I100+"; // no leading simplex count or trailing tail, both meaning zero
        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(100));

        ycTagStr = "20+10A10I20+";

        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(41));

        ycTagStr = "390++"; // empty duplex region and omitted tail, not invalid

        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertTrue(duplexIndelIndices.isEmpty());

        ycTagStr = "++"; // no leading simplex count, empty duplex region, and omitted tail

        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertTrue(duplexIndelIndices.isEmpty());
    }
}
