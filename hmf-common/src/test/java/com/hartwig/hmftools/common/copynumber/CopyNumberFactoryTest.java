package com.hartwig.hmftools.common.copynumber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class CopyNumberFactoryTest {

    @Test
    public void canConvertCNVLine() throws HartwigException {
        final String cnvLineGain = "1 \t 1 \t 2 \t 3 \t gain";
        final String cnvLineLoss = "1 \t 1 \t 2 \t 1 \t loss";

        final CopyNumber cnvGain = CopyNumberFactory.fromCNVLine(cnvLineGain);
        assertEquals("1", cnvGain.chromosome());
        assertEquals(2, cnvGain.start());
        assertEquals(2, cnvGain.end());
        assertEquals(3, cnvGain.value());

        final CopyNumber cnvLoss = CopyNumberFactory.fromCNVLine(cnvLineLoss);
        assertEquals("1", cnvLoss.chromosome());
        assertEquals(2, cnvLoss.start());
        assertEquals(2, cnvLoss.end());
        assertEquals(1, cnvLoss.value());
    }

    @Test(expected = HartwigException.class)
    public void lossGainIsNotRecognized() throws HartwigException {
        final String cnvLine = "1 \t 1 \t 2 \t 3 \t hello world";
        CopyNumberFactory.fromCNVLine(cnvLine);
    }

    @Test(expected = HartwigException.class)
    public void lossGainDoesNotMatchWithValue() throws HartwigException {
        final String cnvLine = "1 \t 1 \t 2 \t 1 \t gain";
        CopyNumberFactory.fromCNVLine(cnvLine);
    }
}