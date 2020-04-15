package com.hartwig.hmftools.patientdb.readers;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideResponseData;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class WidePatientReaderTest {

    @Test
    public void canInterpretPathologySampleId() {
        assertNull(WidePatientReader.extractYearFromPathologySampleId(Strings.EMPTY));
        assertEquals("T1", WidePatientReader.extractYearFromPathologySampleId("T1-00895 P"));

        assertNull(WidePatientReader.extractBiopsyIdFromPathologySampleId(Strings.EMPTY));
        assertNull(WidePatientReader.extractBiopsyIdFromPathologySampleId("weird id"));
        assertEquals("895", WidePatientReader.extractBiopsyIdFromPathologySampleId("T1-00895"));
        assertEquals("895", WidePatientReader.extractBiopsyIdFromPathologySampleId("T1-00895 P"));
        assertEquals("895", WidePatientReader.extractBiopsyIdFromPathologySampleId("T1-00895 P-2"));
    }

    @Test
    public void canDetermineResponse() {
        WideResponseData responseFollowRecist = ImmutableWideResponseData.builder()
                .widePatientId(Strings.EMPTY)
                .timePoint(1)
                .date(LocalDate.parse("2017-01-01"))
                .recistDone(true)
                .recistResponse("PD")
                .noRecistResponse(Strings.EMPTY)
                .noRecistReasonStopTreatment(Strings.EMPTY)
                .noRecistReasonStopTreatmentOther(Strings.EMPTY)
                .build();

        assertEquals(WidePatientReader.determineResponse(responseFollowRecist), "(1) PD");

        WideResponseData responseNotRecist = ImmutableWideResponseData.builder()
                .widePatientId(Strings.EMPTY)
                .timePoint(2)
                .date(LocalDate.parse("2017-01-01"))
                .recistDone(false)
                .recistResponse(Strings.EMPTY)
                .noRecistResponse("stop treatment")
                .noRecistReasonStopTreatment("other, please specify")
                .noRecistReasonStopTreatmentOther(Strings.EMPTY)
                .build();

        assertEquals(WidePatientReader.determineResponse(responseNotRecist), "(2) stop treatment (other, please specify)");
    }

    @Test
    public void canLoadEmptyPatient() {
        WideEcrfModel wideEcrfModel = ImmutableWideEcrfModel.builder()
                .preAvlTreatments(Lists.newArrayList())
                .biopsies(Lists.newArrayList())
                .avlTreatments(Lists.newArrayList())
                .responses(Lists.newArrayList())
                .build();

        WidePatientReader patientReader = new WidePatientReader(wideEcrfModel,
                TestCuratorFactory.tumorLocationCurator(),
                TestCuratorFactory.treatmentCurator());

        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read("ID", "Melanoma", Lists.newArrayList(sample));

        assertNotNull(patient);
    }
}