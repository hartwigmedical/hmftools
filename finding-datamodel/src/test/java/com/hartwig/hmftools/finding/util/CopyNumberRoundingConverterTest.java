package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import junit.framework.TestCase;

public class CopyNumberRoundingConverterTest
{

    private final static double EPSILON = 0.00001;

    public static void assertEqualsDouble(double expected, @Nullable Double actual) {
        assertNotNull(actual);
        assertEquals(expected, actual, EPSILON);
    }

    @Test
    public void testRoundCopyNumber() {
        assertEqualsDouble(2.0, CopyNumberRoundingConverter.roundCopyNumber(2.0136));
        assertEqualsDouble(2.5, CopyNumberRoundingConverter.roundCopyNumber(2.45));
        assertEqualsDouble(2.6, CopyNumberRoundingConverter.roundCopyNumber(2.61));
        assertEqualsDouble(0.0, CopyNumberRoundingConverter.roundCopyNumber(-2.61));
        assertEqualsDouble(114.0, CopyNumberRoundingConverter.roundCopyNumber(113.56));
        assertEqualsDouble(2.4, CopyNumberRoundingConverter.roundCopyNumber(2.44));
        assertNull(CopyNumberRoundingConverter.roundCopyNumberNullable(null));

    }

}