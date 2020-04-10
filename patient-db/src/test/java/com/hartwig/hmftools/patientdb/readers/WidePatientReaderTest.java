package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.*;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideResponseData;
import com.hartwig.hmftools.patientdb.readers.wide.WidePatientReader;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WidePatientReaderTest {
    private static final Logger LOGGER = LogManager.getLogger(WidePatientReaderTest.class);

    @Test
    public void canInterpretDateNL() {
        assertEquals(WidePatientReader.createInterpretDateNL("18-apr-2019").toString(), "2019-04-18");
        assertEquals(WidePatientReader.createInterpretDateNL("17-okt-2018").toString(), "2018-10-17");
    }

    @Test
    public void canInterpretDateEN() {
        assertEquals(WidePatientReader.createInterpretDateEN("21-May-2019").toString(), "2019-05-21");
    }

    @Test
    public void canInterpretDateIC() {
        assertEquals(WidePatientReader.createInterpretDateIC("12/04/2019").toString(), "2019-04-12");
    }

    @Test
    public void canConvertStringToDate() {
        assertEquals(WidePatientReader.convertStringToLocalDate("2018-10-17").toString(), "2018-10-17");
    }

    @Test
    public void canConvertTissueId() {
        assertEquals(WidePatientReader.extractYearOfTissueId("T1-00895 P"), "T1");
        assertEquals(WidePatientReader.extractBiopsyIdOfTissueId("T1-00895 P"), "895");
    }

    @Test
    public void determineResponse() {
        WideResponseData responseFollowRecist = ImmutableWideResponseData.builder()
                .patientId(Strings.EMPTY)
                .timePoint("1")
                .date(Strings.EMPTY)
                .recistNotDone("FALSE")
                .responseAccordingRecist("PD")
                .clinicalDecision(Strings.EMPTY)
                .reasonStopTreatment(Strings.EMPTY)
                .reasonStopTreatmentOther(Strings.EMPTY)
                .build();

        assertEquals(WidePatientReader.determineResponse(responseFollowRecist), "(1) PD");

        WideResponseData responseNotRecist = ImmutableWideResponseData.builder()
                .patientId(Strings.EMPTY)
                .timePoint("2")
                .date(Strings.EMPTY)
                .recistNotDone("WAAR") // TODO change to TRUE
                .responseAccordingRecist(Strings.EMPTY)
                .clinicalDecision("stop treatment")
                .reasonStopTreatment("other, please specify")
                .reasonStopTreatmentOther(Strings.EMPTY)
                .build();

        assertEquals(WidePatientReader.determineResponse(responseNotRecist), "(2) stop treatment (other, please specify)");
    }

    @Test
    public void convertGender() {
        assertEquals(WidePatientReader.convertGender("1"), "Male");
        assertEquals(WidePatientReader.convertGender("2"), "Female");
        assertEquals(WidePatientReader.convertGender(""), Strings.EMPTY);
    }

    //    @Test
    //    @Ignore
    //    public void canLoadEmptyPatient() {
    //
    //
    //        WideEcrfModel wideEcrfModel;
    //
    //        wideEcrfModel = ImmutableWideEcrfModel.builder()
    //                .preTreatments(Lists.newArrayList())
    //                .biopsies(Lists.newArrayList())
    //                .treatments(Lists.newArrayList())
    //                .responses(Lists.newArrayList())
    //                .build();
    //
    //        WidePatientReader patientReader = new WidePatientReader(wideEcrfModel, TestCuratorFactory.tumorLocationCurator(),
    //                TestCuratorFactory.biopsySiteCurator(),
    //                TestCuratorFactory.treatmentCurator());
    //
    //        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();
    //
    //        Patient patient = patientReader.read("ID", "Melanoma", Lists.newArrayList(sample));
    //
    //        assertNotNull(patient);
    //    }
}