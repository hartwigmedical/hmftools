package com.hartwig.hmftools.ckb.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class DateConverterTest {

    @Test
    public void canConvertString() {
        String date = "03/29/2017";
        assertEquals(date, DateConverter.FORMAT.format(DateConverter.toDate(date)));

        assertNull(DateConverter.toDate(null));
    }

    @Test (expected = IllegalStateException.class)
    public void wrongDateLeadsToException() {
        DateConverter.toDate("not a date");
    }
}