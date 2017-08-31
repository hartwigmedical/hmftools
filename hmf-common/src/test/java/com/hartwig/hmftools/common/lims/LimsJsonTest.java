package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class LimsJsonTest {

    @Test
    public void canLoadJsonModel() throws IOException, EmptyFileException {
        final String limsJson = Resources.getResource("lims").getPath() + File.separator + "lims.json";
        final LimsJsonModel jsonLims = LimsJsonModel.readModelFromFile(limsJson);
        assertEquals(LocalDate.parse("2016-01-02"), jsonLims.samplingDateForSample("SAMP01010003R"));
        assertEquals(LocalDate.parse("2016-01-03"), jsonLims.arrivalDateForSample("SAMP01010003R"));
        assertNull(jsonLims.tumorPercentageForSample("SAMP01010003R"));
        assertEquals("CSB000000", jsonLims.bloodBarcodeForSample("SAMP01010003R"));
        assertEquals("CSB000000", jsonLims.barcodeForSample("SAMP01010003R"));
        assertEquals(LocalDate.parse("2016-01-03"), jsonLims.bloodArrivalDateForSample("SAMP01010003R"));

        assertEquals("CSB000000", jsonLims.bloodBarcodeForSample("SAMP01010003T"));
        assertEquals("FC0000001", jsonLims.barcodeForSample("SAMP01010003T"));
        final Double tumorPercentage = jsonLims.tumorPercentageForSample("SAMP01010003T");
        assertNotNull(tumorPercentage);
        assertEquals(.3, tumorPercentage, 0.001);
        assertEquals(LocalDate.parse("2016-02-05"), jsonLims.arrivalDateForSample("SAMP01010003T"));
        assertEquals(LocalDate.parse("2016-01-04"), jsonLims.samplingDateForSample("SAMP01010003T"));
        assertEquals(LocalDate.parse("2016-01-03"), jsonLims.bloodArrivalDateForSample("SAMP01010003T"));
    }
}
