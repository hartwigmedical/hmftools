package com.hartwig.hmftools.ckb.util;

import static org.junit.Assert.*;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Locale;

import org.junit.Test;

public class DateConverterTest {

    @Test
    public void canConvertString() throws ParseException {
        String date = "10/02/2020";
        SimpleDateFormat format = new SimpleDateFormat("dd/MM/yyyy", Locale.ENGLISH);
        assertEquals(date, format.format(DateConverter.convertDate(date)));
    }

}