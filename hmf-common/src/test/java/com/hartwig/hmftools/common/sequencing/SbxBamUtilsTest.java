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
        String ycTagStr = "0+100+0";
        List<Integer> duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertTrue(duplexIndelIndices.isEmpty());

        ycTagStr = "10+100I100+0";
        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(110));

        ycTagStr = "+100I100+0"; // no leading simplex count, meaning zero
        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(100));

        ycTagStr = "20+10A10I20+";

        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(41));
    }
}
