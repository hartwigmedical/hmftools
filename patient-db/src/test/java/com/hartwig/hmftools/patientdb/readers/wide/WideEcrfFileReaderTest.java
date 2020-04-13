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
    public void canReadWidePreAvlTreatments() throws IOException {
        List<WidePreAVLTreatmentData> preTreatments =
                WideEcrfFileReader.readPreAvlTreatmentData(WIDE_TEST_DIR + File.separator + "wide_pre_avl_treatments.csv");

        assertEquals(2, preTreatments.size());
    }

    @Test
    public void canReadWideBiopsies() throws IOException {
        List<WideBiopsyData> biopsies =
                WideEcrfFileReader.readBiopsyData(WIDE_TEST_DIR + File.separator + "wide_biopsies.csv");

        assertEquals(1, biopsies.size());
    }

    @Test
    public void canReadWideAvlTreatments() throws IOException {
        List<WideAvlTreatmentData> treatments =
                WideEcrfFileReader.readAvlTreatmentData(WIDE_TEST_DIR + File.separator + "wide_avl_treatments.csv");

        assertEquals(2, treatments.size());
    }

    @Test
    public void canReadWideResponses() throws IOException {
        List<WideResponseData> responses =
                WideEcrfFileReader.readResponseData(WIDE_TEST_DIR + File.separator + "wide_responses.csv");

        assertEquals(2, responses.size());
    }

    @Test
    public void canInterpretDateNL() {
        assertEquals(WideEcrfFileReader.interpretDateNL("18-apr-2019").toString(), "2019-04-18");
        assertEquals(WideEcrfFileReader.interpretDateNL("17-okt-2018").toString(), "2018-10-17");
    }

    @Test
    public void canInterpretDateEN() {
        assertEquals(WideEcrfFileReader.interpretDateEN("21-May-2019").toString(), "2019-05-21");
    }
}