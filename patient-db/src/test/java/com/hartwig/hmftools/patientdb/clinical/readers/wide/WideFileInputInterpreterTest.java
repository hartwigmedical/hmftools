package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;

import org.junit.Test;

public class WideFileInputInterpreterTest
{

    @Test
    public void canInterpretDateIC()
    {
        assertEquals(LocalDate.parse("2019-04-12"), WideFileInputInterpreter.interpretDateIC("12/04/2019"));
    }

    @Test
    public void canInterpretDateEN()
    {
        assertEquals(LocalDate.parse("2019-05-21"), WideFileInputInterpreter.interpretDateEN("21-May-2019"));
    }

    @Test
    public void canInterpretBirthyear()
    {
        assertEquals(1960, (int)WideFileInputInterpreter.interpretBirthYear("1960"));
        assertNull(WideFileInputInterpreter.interpretBirthYear(""));
    }

    @Test
    public void canConvertGender()
    {
        assertEquals("male", WideFileInputInterpreter.convertGender("1"));
        assertEquals("female", WideFileInputInterpreter.convertGender("2"));
        assertNull(WideFileInputInterpreter.convertGender(""));
    }

    @Test
    public void canConvertParticipatesInOtherTrials()
    {
        assertTrue(WideFileInputInterpreter.convertParticipatesInOtherTrials("Y"));
        assertFalse(WideFileInputInterpreter.convertParticipatesInOtherTrials("N"));
        assertNull(WideFileInputInterpreter.convertParticipatesInOtherTrials(""));
    }

    @Test
    public void canConvertHasSuccessfulReport()
    {
        assertTrue(WideFileInputInterpreter.convertHasReceivedSuccessfulReport("yes"));
        assertFalse(WideFileInputInterpreter.convertHasReceivedSuccessfulReport("no"));
        assertNull(WideFileInputInterpreter.convertHasReceivedSuccessfulReport(""));
    }

    @Test
    public void canConvertWideID()
    {
        assertEquals("WIDE01018888", WideFileInputInterpreter.toWideID("WIDE-01-01-8888"));
    }
}