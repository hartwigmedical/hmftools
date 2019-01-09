package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDate;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canBuildLimsFromTestData() throws IOException {
        assertNotNull(LimsFactory.fromLimsDirectory(LIMS_DIRECTORY));
    }

    @Test(expected = IOException.class)
    public void exceptionWhenJsonFileDoesNotExist() throws IOException {
        LimsFactory.fromLimsDirectory("Does not exist");
    }

    @Test
    public void readCorrectlyFromJsonFile() throws FileNotFoundException {
        final Map<String, LimsJsonData> dataPerSample = LimsFactory.readLimsJson(LIMS_DIRECTORY + File.separator + "lims.json");

        final LimsJsonData refData = dataPerSample.get("SAMP01010003R");
        assertEquals("2016-01-02", refData.samplingDateString());
        assertEquals("2016-01-03", refData.arrivalDateString());
        assertEquals("NA", refData.tumorPercentageString());
        assertEquals("PREP013V23-QC037V20-SEQ008V25", refData.labProcedures());

        final LimsJsonData tumorData = dataPerSample.get("SAMP01010003T");
        assertEquals("2016-01-04", tumorData.samplingDateString());
        assertEquals("2016-02-05", tumorData.arrivalDateString());
        assertEquals("30", tumorData.tumorPercentageString());
        assertEquals("N/A", tumorData.labProcedures());
    }

    @Test
    public void readCorrectlyFromPreLIMSArrivalDateFile() throws IOException {
        final Map<String, LocalDate> preLIMSArrivalDates =
                LimsFactory.readPreLIMSArrivalDateCsv(LIMS_DIRECTORY + File.separator + "pre_lims_arrival_dates.csv");

        assertNull(preLIMSArrivalDates.get("SAMP01010001T"));
        assertEquals(LimsTestUtil.toDate("2017-01-01"), preLIMSArrivalDates.get("SAMP01010003R"));
        assertEquals(LimsTestUtil.toDate("2017-01-05"), preLIMSArrivalDates.get("SAMP01010004T"));
        assertNull(preLIMSArrivalDates.get("SAMP01010005T"));
        assertNull(preLIMSArrivalDates.get("SAMP01010006R"));
        assertNull(preLIMSArrivalDates.get("DoesNotExist"));
    }
}