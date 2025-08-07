package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndels;

import static org.junit.Assert.assertEquals;
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
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(100, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsNoDuplexIndelsWithSimplexHead()
    {
        String ycTagStr = "10-100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(110, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsSingleDuplexIndel()
    {
        String ycTagStr = "0-100I100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(201, duplexIndels.size());
        assertTrue(duplexIndels.subList(0, 100).stream().noneMatch(x -> x));
        assertTrue(duplexIndels.get(100));
        assertTrue(duplexIndels.subList(101, 201).stream().noneMatch(x -> x));
    }

    @Test
    public void testGetDuplexIndelsSingleNonDuplexIndel()
    {
        String ycTagStr = "0-100A100-0";
        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);

        assertEquals(201, duplexIndels.size());
        assertTrue(duplexIndels.stream().noneMatch(x -> x));
    }
}
