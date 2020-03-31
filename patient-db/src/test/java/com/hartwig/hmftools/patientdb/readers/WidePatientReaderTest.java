package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.*;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.Locale;

import org.junit.Test;

public class WidePatientReaderTest {

        @Test
        public void canInterpretDate() {

            DateTimeFormatter inputFormatter = new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH);
            LocalDate localDate = LocalDate.parse("18-apr-2019", inputFormatter);

            DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
            String formattedString = localDate.format(outputFormatter);
            LocalDate date = LocalDate.parse(formattedString);

            assertEquals(date.toString(), "2019-04-18");
        }
}