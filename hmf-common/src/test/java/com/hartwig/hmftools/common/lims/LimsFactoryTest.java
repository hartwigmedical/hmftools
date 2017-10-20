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

    private static final String LIMS_RESOURCES = Resources.getResource("lims").getPath();
    private static final String LIMS_JSON = LIMS_RESOURCES + File.separator + "lims_test.json";
    private static final String PRE_HMF_ARRIVAL_DATES_CSV = LIMS_RESOURCES + File.separator + "pre_hmf_arrival_dates.csv";

    @Test
    public void canBuildLimsFromTestData() throws IOException {
        assertNotNull(LimsFactory.fromLimsJson(LIMS_JSON));
        assertNotNull(LimsFactory.fromLimsJsonWithPreHMFArrivalDates(LIMS_JSON, PRE_HMF_ARRIVAL_DATES_CSV));
    }

    @Test(expected = FileNotFoundException.class)
    public void exceptionWhenJsonFileDoesNotExist() throws FileNotFoundException {
        LimsFactory.fromLimsJson("Does not exist");
    }

    @Test
    public void readCorrectlyFromJsonFile() throws FileNotFoundException {
        final Map<String, LimsJsonData> dataPerSample = LimsFactory.readLimsJson(LIMS_JSON);

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
    public void readCorrectlyFromPreHMFArrivalDateFile() throws IOException {
        final Map<String, LocalDate> preHMFArrivalDates = LimsFactory.readPreHMFArrivalDateCsv(PRE_HMF_ARRIVAL_DATES_CSV);

        assertNull(preHMFArrivalDates.get("SAMP01010001T"));
        assertEquals(LimsTestUtil.toDate("2017-01-01"), preHMFArrivalDates.get("SAMP01010003R"));
        assertEquals(LimsTestUtil.toDate("2017-01-05"), preHMFArrivalDates.get("SAMP01010004T"));
        assertNull(preHMFArrivalDates.get("SAMP01010005T"));
        assertNull(preHMFArrivalDates.get("SAMP01010006R"));
        assertNull(preHMFArrivalDates.get("DoesNotExist"));
    }
}