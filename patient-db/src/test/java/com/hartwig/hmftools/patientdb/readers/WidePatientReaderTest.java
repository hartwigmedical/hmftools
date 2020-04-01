package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.*;

import org.junit.Test;

public class WidePatientReaderTest {

        @Test
        public void canInterpretDate() {

            assertEquals(WidePatientReader.createInterpretDate("18-apr-2019").toString(), "2019-04-18");
            assertEquals(WidePatientReader.createInterpretDate("17-okt-2018").toString(), "2018-10-17");

        }
}