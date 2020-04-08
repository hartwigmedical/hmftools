package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.*;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import com.hartwig.hmftools.patientdb.readers.wide.WidePatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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