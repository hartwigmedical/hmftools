package com.hartwig.hmftools.ckb.util;

import static org.junit.Assert.assertEquals;

import java.time.format.DateTimeFormatter;
import java.util.Locale;

import org.junit.Test;

public class DateConverterTest {

    @Test
    public void canConvertString() {
        String date = "10/02/2020";
        DateTimeFormatter format = DateTimeFormatter.ofPattern("dd/MM/yyyy", Locale.ENGLISH);
        assertEquals(date, format.format(DateConverter.toDate(date)));
    }
}