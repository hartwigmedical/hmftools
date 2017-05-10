package com.hartwig.hmftools.common.freec;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.exception.HartwigException;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class FreecCopyNumberFactoryTest {

    private static final String BASE_PATH = Resources.getResource("copynumber").getPath();
    private static final String SAMPLE = "sample";

    @Test
    public void canLoadCNVFile() throws IOException, HartwigException {
        final List<CopyNumber> copyNumbers = FreecCopyNumberFactory.loadCNV(BASE_PATH, SAMPLE);
        assertEquals(2, copyNumbers.size());
    }

    @Test
    public void canConvertCNVLine() throws HartwigException {
        final String cnvLineGain = "1 \t 1 \t 2 \t 3 \t gain \t AAB \t 27.30 \t germline";
        final String cnvLineLoss = "1 \t 1 \t 2 \t 1 \t loss";

        final FreecCopyNumber cnvGain = FreecCopyNumberFactory.fromCNVLine(cnvLineGain);
        assertEquals("1", cnvGain.chromosome());
        assertEquals(2, cnvGain.start());
        assertEquals(2, cnvGain.end());
        assertEquals(3, cnvGain.value());
        assertEquals("AAB", cnvGain.genotype());
        assertEquals(FreecCopyNumber.Status.GERMLINE, cnvGain.status());

        final FreecCopyNumber cnvLoss = FreecCopyNumberFactory.fromCNVLine(cnvLineLoss);
        assertEquals("1", cnvLoss.chromosome());
        assertEquals(2, cnvLoss.start());
        assertEquals(2, cnvLoss.end());
        assertEquals(1, cnvLoss.value());
        assertEquals("-", cnvLoss.genotype());
        assertEquals(FreecCopyNumber.Status.UNKNOWN, cnvLoss.status());
    }

    @Test(expected = HartwigException.class)
    public void lossGainIsNotRecognized() throws HartwigException {
        final String cnvLine = "1 \t 1 \t 2 \t 3 \t hello world";
        FreecCopyNumberFactory.fromCNVLine(cnvLine);
    }

    @Test(expected = HartwigException.class)
    public void lossGainDoesNotMatchWithValue() throws HartwigException {
        final String cnvLine = "1 \t 1 \t 2 \t 1 \t gain";
        FreecCopyNumberFactory.fromCNVLine(cnvLine);
    }
}