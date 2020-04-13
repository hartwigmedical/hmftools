package com.hartwig.hmftools.patientdb.readers.wide;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class WideEcrfFileReaderTest {

    private static final String WIDE_TEST_DIR = Resources.getResource("wide").getPath();

    @Test
    public void canReadWidePreAvlTreatments() throws IOException {
        List<WidePreAvlTreatmentData> preTreatments =
                WideEcrfFileReader.readPreAvlTreatments(WIDE_TEST_DIR + File.separator + "wide_pre_avl_treatments.csv");

        assertEquals(2, preTreatments.size());
    }

    @Test
    public void canReadWideBiopsies() throws IOException {
        List<WideBiopsyData> biopsies = WideEcrfFileReader.readBiopsies(WIDE_TEST_DIR + File.separator + "wide_biopsies.csv");

        assertEquals(1, biopsies.size());
    }

    @Test
    public void canReadWideAvlTreatments() throws IOException {
        List<WideAvlTreatmentData> treatments =
                WideEcrfFileReader.readAvlTreatments(WIDE_TEST_DIR + File.separator + "wide_avl_treatments.csv");

        assertEquals(2, treatments.size());
    }

    @Test
    public void canReadWideResponses() throws IOException {
        List<WideResponseData> responses = WideEcrfFileReader.readResponses(WIDE_TEST_DIR + File.separator + "wide_responses.csv");

        assertEquals(2, responses.size());
    }

    @Test
    public void canReadWideFiveDays() throws IOException {
        List<WideFiveDays> fiveDays = WideEcrfFileReader.readFiveDays(WIDE_TEST_DIR + File.separator + "wide_five_days.csv");

        assertEquals(3, fiveDays.size());

        assertEquals("WIDE01018888", fiveDays.get(0).patientId());
        assertTrue(fiveDays.get(0).dataIsAvailable());
        assertEquals(LocalDate.parse("2018-03-12"), fiveDays.get(0).informedConsentDate());
        assertEquals("female", fiveDays.get(0).gender());
        assertEquals(1940, (int) fiveDays.get(0).birthYear());
        assertEquals(LocalDate.parse("2018-04-11"), fiveDays.get(0).biopsyDate());
        assertEquals("peritoneum", fiveDays.get(0).biopsySite());
        assertEquals("peritoneum", fiveDays.get(0).sampleTissue());
        assertEquals("biopt", fiveDays.get(0).sampleType());
        assertEquals("N18WGS-2", fiveDays.get(0).studyCode());
        assertFalse(fiveDays.get(0).participatesInOtherTrials());
        assertTrue(fiveDays.get(0).otherTrialCodes().isEmpty());
        assertTrue(fiveDays.get(0).otherTrialStartDates().isEmpty());

        assertEquals("WIDE01018888", fiveDays.get(1).patientId());
        assertTrue(fiveDays.get(1).dataIsAvailable());
        assertEquals(LocalDate.parse("2018-03-12"), fiveDays.get(1).informedConsentDate());
        assertEquals("female", fiveDays.get(1).gender());
        assertEquals(1940, (int) fiveDays.get(1).birthYear());
        assertEquals(LocalDate.parse("2018-05-12"), fiveDays.get(1).biopsyDate());
        assertEquals("lever", fiveDays.get(1).biopsySite());
        assertEquals("lever neoplasie", fiveDays.get(1).sampleTissue());
        assertEquals("biopt", fiveDays.get(1).sampleType());
        assertEquals("N18WGS-2", fiveDays.get(1).studyCode());
        assertFalse(fiveDays.get(1).participatesInOtherTrials());
        assertTrue(fiveDays.get(1).otherTrialCodes().isEmpty());
        assertTrue(fiveDays.get(1).otherTrialStartDates().isEmpty());

        assertEquals("WIDE01019999", fiveDays.get(2).patientId());
        assertFalse(fiveDays.get(2).dataIsAvailable());
        assertNull(fiveDays.get(2).informedConsentDate());
        assertNull(fiveDays.get(2).gender());
        assertNull(fiveDays.get(2).birthYear());
        assertNull(fiveDays.get(2).biopsyDate());
        assertTrue(fiveDays.get(2).biopsySite().isEmpty());
        assertTrue(fiveDays.get(2).sampleTissue().isEmpty());
        assertTrue(fiveDays.get(2).sampleType().isEmpty());
        assertTrue(fiveDays.get(2).studyCode().isEmpty());
        assertNull(fiveDays.get(2).participatesInOtherTrials());
        assertTrue(fiveDays.get(2).otherTrialCodes().isEmpty());
        assertTrue(fiveDays.get(2).otherTrialStartDates().isEmpty());
    }

    @Test
    public void canInterpretDateIC() {
        assertEquals(LocalDate.parse("2019-04-12"), WideEcrfFileReader.interpretDateIC("12/04/2019"));
    }

    @Test
    public void canInterpretDateNL() {
        assertEquals(LocalDate.parse("2019-04-18"), WideEcrfFileReader.interpretDateNL("18-apr-2019"));
        assertEquals(LocalDate.parse("2018-10-17"), WideEcrfFileReader.interpretDateNL("17-okt-2018"));
    }

    @Test
    public void canInterpretDateEN() {
        assertEquals(LocalDate.parse("2019-05-21"), WideEcrfFileReader.interpretDateEN("21-May-2019"));
    }

    @Test
    public void canConvertGender() {
        assertEquals("male", WideEcrfFileReader.convertGender("1"));
        assertEquals("female", WideEcrfFileReader.convertGender("2"));
        assertNull(WideEcrfFileReader.convertGender(""));
    }

    @Test
    public void canConvertParticipatesInOtherTrials() {
        assertTrue("male", WideEcrfFileReader.convertParticipatesInOtherTrials("Y"));
        assertFalse("female", WideEcrfFileReader.convertParticipatesInOtherTrials("N"));
        assertNull(WideEcrfFileReader.convertParticipatesInOtherTrials(""));
    }
}