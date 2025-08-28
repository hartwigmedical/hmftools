package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndelIndices;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndels;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.List;

import org.junit.Test;

public class SbxBamUtilsTest
{
    @Test
    public void testGetDuplexIndelsNoDuplexIndels()
    {
        String ycTagStr = "0-100-0";
        List<Integer> duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertNull(duplexIndelIndices);

        ycTagStr = "10-100I100-0";
        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(110));

        ycTagStr = "20-10A10I20-0";

        duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        assertEquals(1, duplexIndelIndices.size());
        assertTrue(duplexIndelIndices.contains(41));
    }
}
