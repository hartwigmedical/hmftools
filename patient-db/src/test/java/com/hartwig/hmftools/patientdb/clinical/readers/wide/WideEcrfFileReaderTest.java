package com.hartwig.hmftools.patientdb.clinical.readers.wide;

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

public class WideEcrfFileReaderTest
{

    private static final String WIDE_TEST_DIR = Resources.getResource("wide").getPath();

    @Test
    public void canReadWidePreAvlTreatments() throws IOException
    {
        List<WidePreAvlTreatmentData> preTreatments =
                WideEcrfFileReader.readPreAvlTreatments(WIDE_TEST_DIR + File.separator + "wide_pre_avl_treatments.csv");

        assertEquals(3, preTreatments.size());

        assertEquals("WIDE01018888", preTreatments.get(0).widePatientId());
        assertFalse(preTreatments.get(0).hasPreviousTherapy());
        assertTrue(preTreatments.get(0).drug1().isEmpty());
        assertTrue(preTreatments.get(0).drug2().isEmpty());
        assertTrue(preTreatments.get(0).drug3().isEmpty());
        assertTrue(preTreatments.get(0).drug4().isEmpty());
        assertNull(preTreatments.get(0).lastSystemicTherapyDate());

        assertEquals("WIDE01019999", preTreatments.get(1).widePatientId());
        assertTrue(preTreatments.get(1).hasPreviousTherapy());
        assertEquals("erlotinib", preTreatments.get(1).drug1());
        assertTrue(preTreatments.get(1).drug2().isEmpty());
        assertTrue(preTreatments.get(1).drug3().isEmpty());
        assertTrue(preTreatments.get(1).drug4().isEmpty());
        assertEquals(LocalDate.parse("2018-12-17"), preTreatments.get(1).lastSystemicTherapyDate());

        assertEquals("WIDE01019999", preTreatments.get(2).widePatientId());
        assertTrue(preTreatments.get(2).hasPreviousTherapy());
        assertEquals("carboplatin,others", preTreatments.get(2).drug1());
        assertEquals("pemetrexed", preTreatments.get(2).drug2());
        assertEquals("irinotecan", preTreatments.get(2).drug3());
        assertEquals("bevacizumab", preTreatments.get(2).drug4());
        assertEquals(LocalDate.parse("2018-12-17"), preTreatments.get(2).lastSystemicTherapyDate());
    }

    @Test
    public void canReadWideBiopsies() throws IOException
    {
        List<WideBiopsyData> biopsies = WideEcrfFileReader.readBiopsies(WIDE_TEST_DIR + File.separator + "wide_biopsies.csv");

        assertEquals(2, biopsies.size());

        assertEquals("WIDE01018888", biopsies.get(0).widePatientId());
        assertEquals("T18-00001", biopsies.get(0).pathologySampleId());
        assertEquals(LocalDate.parse("2018-01-18"), biopsies.get(0).biopsyDate());
        assertNull(biopsies.get(0).hasReceivedSuccessfulReport());

        assertEquals("WIDE01019999", biopsies.get(1).widePatientId());
        assertEquals("T19-00001", biopsies.get(1).pathologySampleId());
        assertEquals(LocalDate.parse("2019-04-18"), biopsies.get(1).biopsyDate());
        assertTrue(biopsies.get(1).hasReceivedSuccessfulReport());
    }

    @Test
    public void canReadWideAvlTreatments() throws IOException
    {
        List<WideAvlTreatmentData> treatments =
                WideEcrfFileReader.readAvlTreatments(WIDE_TEST_DIR + File.separator + "wide_avl_treatments.csv");

        assertEquals(4, treatments.size());

        assertEquals("WIDE01018888", treatments.get(0).widePatientId());
        assertEquals("L01BA04", treatments.get(0).drugCode());
        assertEquals("pemetrexed", treatments.get(0).drug());
        assertEquals(LocalDate.parse("2018-03-03"), treatments.get(0).startDate());
        assertEquals(LocalDate.parse("2018-07-20"), treatments.get(0).endDate());

        assertEquals("WIDE01019999", treatments.get(1).widePatientId());
        assertEquals("L01XE03", treatments.get(1).drugCode());
        assertEquals("erlotinib", treatments.get(1).drug());
        assertEquals(LocalDate.parse("2018-10-06"), treatments.get(1).startDate());
        assertNull(treatments.get(1).endDate());

        assertEquals("WIDE01019999", treatments.get(2).widePatientId());
        assertEquals("L01XE03", treatments.get(2).drugCode());
        assertEquals("erlotinib", treatments.get(2).drug());
        assertEquals(LocalDate.parse("2018-08-18"), treatments.get(2).startDate());
        assertNull(treatments.get(2).endDate());
    }

    @Test
    public void canReadWideAvlTreatmentsSemi() throws IOException
    {
        List<WideAvlTreatmentData> treatments =
                WideEcrfFileReader.readAvlTreatments(WIDE_TEST_DIR + File.separator + "wide_avl_treatments_semi.csv");

        assertEquals(4, treatments.size());

        assertEquals("WIDE01018888", treatments.get(0).widePatientId());
        assertEquals("L01BA04", treatments.get(0).drugCode());
        assertEquals("pemetrexed", treatments.get(0).drug());
        assertEquals(LocalDate.parse("2018-03-03"), treatments.get(0).startDate());
        assertEquals(LocalDate.parse("2018-07-20"), treatments.get(0).endDate());

        assertEquals("WIDE01019999", treatments.get(1).widePatientId());
        assertEquals("L01XE03", treatments.get(1).drugCode());
        assertEquals("erlotinib", treatments.get(1).drug());
        assertEquals(LocalDate.parse("2018-10-06"), treatments.get(1).startDate());
        assertNull(treatments.get(1).endDate());

        assertEquals("WIDE01019999", treatments.get(2).widePatientId());
        assertEquals("L01XE03", treatments.get(2).drugCode());
        assertEquals("erlotinib", treatments.get(2).drug());
        assertEquals(LocalDate.parse("2018-09-18"), treatments.get(2).startDate());
        assertNull(treatments.get(2).endDate());
    }

    @Test
    public void canReadWideResponses() throws IOException
    {
        List<WideResponseData> responses = WideEcrfFileReader.readResponses(WIDE_TEST_DIR + File.separator + "wide_responses.csv");

        assertEquals(3, responses.size());

        assertEquals("WIDE01018888", responses.get(0).widePatientId());
        assertEquals(1, responses.get(0).timePoint());
        assertEquals(LocalDate.parse("2019-05-05"), responses.get(0).date());
        assertFalse(responses.get(0).recistDone());
        assertTrue(responses.get(0).recistResponse().isEmpty());
        assertTrue(responses.get(0).noRecistResponse().isEmpty());
        assertTrue(responses.get(0).noRecistReasonStopTreatment().isEmpty());
        assertTrue(responses.get(0).noRecistReasonStopTreatmentOther().isEmpty());

        assertEquals("WIDE01018888", responses.get(1).widePatientId());
        assertEquals(2, responses.get(1).timePoint());
        assertEquals(LocalDate.parse("2019-06-28"), responses.get(1).date());
        assertFalse(responses.get(1).recistDone());
        assertTrue(responses.get(1).recistResponse().isEmpty());
        assertEquals("stop treatment", responses.get(1).noRecistResponse());
        assertEquals("radiological progression of non-RECIST lesions", responses.get(1).noRecistReasonStopTreatment());
        assertTrue(responses.get(1).noRecistReasonStopTreatmentOther().isEmpty());

        assertEquals("WIDE01019999", responses.get(2).widePatientId());
        assertEquals(1, responses.get(2).timePoint());
        assertEquals(LocalDate.parse("2019-09-24"), responses.get(2).date());
        assertTrue(responses.get(2).recistDone());
        assertEquals("PD", responses.get(2).recistResponse());
        assertNull(responses.get(2).noRecistResponse());
        assertNull(responses.get(2).noRecistReasonStopTreatment());
        assertNull(responses.get(2).noRecistReasonStopTreatmentOther());
    }

    @Test
    public void canReadWideFiveDays() throws IOException
    {
        List<WideFiveDays> fiveDays = WideEcrfFileReader.readFiveDays(WIDE_TEST_DIR + File.separator + "wide_five_days.csv");

        assertEquals(3, fiveDays.size());

        assertEquals("WIDE01018888", fiveDays.get(0).widePatientId());
        assertTrue(fiveDays.get(0).dataIsAvailable());
        assertEquals(LocalDate.parse("2018-03-12"), fiveDays.get(0).informedConsentDate());
        assertEquals("female", fiveDays.get(0).gender());
        assertEquals(1940, (int)fiveDays.get(0).birthYear());
        assertEquals(LocalDate.parse("2018-04-11"), fiveDays.get(0).biopsyDate());
        assertEquals("peritoneum", fiveDays.get(0).biopsySite());
        assertEquals("peritoneum", fiveDays.get(0).sampleTissue());
        assertEquals("biopt", fiveDays.get(0).sampleType());
        assertEquals("N18WGS-2", fiveDays.get(0).studyCode());
        assertFalse(fiveDays.get(0).participatesInOtherTrials());
        assertTrue(fiveDays.get(0).otherTrialCodes().isEmpty());
        assertTrue(fiveDays.get(0).otherTrialStartDates().isEmpty());

        assertEquals("WIDE01018888", fiveDays.get(1).widePatientId());
        assertTrue(fiveDays.get(1).dataIsAvailable());
        assertEquals(LocalDate.parse("2018-03-12"), fiveDays.get(1).informedConsentDate());
        assertEquals("female", fiveDays.get(1).gender());
        assertEquals(1940, (int)fiveDays.get(1).birthYear());
        assertEquals(LocalDate.parse("2018-05-12"), fiveDays.get(1).biopsyDate());
        assertEquals("lever", fiveDays.get(1).biopsySite());
        assertEquals("lever neoplasie", fiveDays.get(1).sampleTissue());
        assertEquals("biopt", fiveDays.get(1).sampleType());
        assertEquals("N18WGS-2", fiveDays.get(1).studyCode());
        assertFalse(fiveDays.get(1).participatesInOtherTrials());
        assertTrue(fiveDays.get(1).otherTrialCodes().isEmpty());
        assertTrue(fiveDays.get(1).otherTrialStartDates().isEmpty());

        assertEquals("WIDE01019999", fiveDays.get(2).widePatientId());
        assertFalse(fiveDays.get(2).dataIsAvailable());
        assertNull(fiveDays.get(2).informedConsentDate());
        assertNull(fiveDays.get(2).gender());
        assertNull(fiveDays.get(2).birthYear());
        assertNull(fiveDays.get(2).biopsyDate());
        assertNull(fiveDays.get(2).biopsySite());
        assertNull(fiveDays.get(2).sampleTissue());
        assertNull(fiveDays.get(2).sampleType());
        assertNull(fiveDays.get(2).studyCode());
        assertNull(fiveDays.get(2).participatesInOtherTrials());
        assertNull(fiveDays.get(2).otherTrialCodes());
        assertNull(fiveDays.get(2).otherTrialStartDates());
    }
}