package com.hartwig.hmftools.pavereverse.gene;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Test;

public class AnnotatedExonTest extends ReversePaveTestBase
{
    private final ExonData exonData = new ExonData(123, 105, 115, 2, 0, 0);
    private final AnnotatedExon annotatedExon = new AnnotatedExon(10, 20, exonData);
    private final AnnotatedExon annotatedExonRS = new AnnotatedExon(10, 20, exonData, false);

    @Test
    public void containsTest()
    {
        assertTrue(annotatedExon.contains(10));
        assertTrue(annotatedExonRS.contains(10));
        assertTrue(annotatedExon.contains(15));
        assertTrue(annotatedExonRS.contains(15));
        assertTrue(annotatedExon.contains(20));
        assertTrue(annotatedExonRS.contains(20));
        assertFalse(annotatedExon.contains(9));
        assertFalse(annotatedExonRS.contains(9));
        assertFalse(annotatedExon.contains(21));
        assertFalse(annotatedExonRS.contains(21));
    }

    @Test
    public void getAbsolutePosition()
    {
        assertEquals(105, annotatedExon.getAbsolutePosition(10));
        assertEquals(106, annotatedExon.getAbsolutePosition(11));
        assertEquals(115, annotatedExon.getAbsolutePosition(20));
    }

    @Test
    public void getAbsolutePositionReverseStrand()
    {
        assertEquals(115, annotatedExonRS.getAbsolutePosition(10));
        assertEquals(114, annotatedExonRS.getAbsolutePosition(11));
        assertEquals(110, annotatedExonRS.getAbsolutePosition(15));
        assertEquals(106, annotatedExonRS.getAbsolutePosition(19));
        assertEquals(105, annotatedExonRS.getAbsolutePosition(20));
    }
}
