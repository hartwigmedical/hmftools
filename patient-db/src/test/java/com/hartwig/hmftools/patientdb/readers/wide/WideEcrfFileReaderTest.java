package com.hartwig.hmftools.patientdb.readers.wide;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class WideEcrfFileReaderTest {

    private static final String WIDE_TEST_DIR = Resources.getResource("wide").getPath();

    @Test
    public void canReadWidePreTreatments() throws IOException {
        List<WidePreTreatmentData> preTreatments =
                WideEcrfFileReader.readPreTreatmentData(WIDE_TEST_DIR + File.separator + "wide_pre_treatments.csv");

        assertEquals(2, preTreatments.size());
    }

    @Test
    public void canReadWideBiopsies() throws IOException {
        List<WideBiopsyData> biopsies =
                WideEcrfFileReader.readBiopsyData(WIDE_TEST_DIR + File.separator + "wide_biopsies.csv");

        assertEquals(1, biopsies.size());
    }

    @Test
    public void canReadWideTreatments() throws IOException {
        List<WideTreatmentData> treatments =
                WideEcrfFileReader.readTreatmentData(WIDE_TEST_DIR + File.separator + "wide_treatments.csv");

        assertEquals(2, treatments.size());
    }

    @Test
    public void canReadPreResponses() throws IOException {
        List<WideResponseData> responses =
                WideEcrfFileReader.readResponseData(WIDE_TEST_DIR + File.separator + "wide_responses.csv");

        assertEquals(2, responses.size());
    }
}