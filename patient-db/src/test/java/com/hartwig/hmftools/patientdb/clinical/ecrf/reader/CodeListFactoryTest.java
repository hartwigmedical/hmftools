package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CodeListFactoryTest {

    @Test
    public void canExtractValuesFromStrings() {
        assertEquals("x", CodeListFactory.fromText("1=x"));
        assertEquals("hi", CodeListFactory.fromText("hi"));
    }
}