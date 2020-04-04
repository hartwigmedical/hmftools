package com.hartwig.hmftools.patientdb.readers;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.*;

import java.time.LocalDate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WidePatientReader;

import org.junit.Test;

public class WidePatientReaderTest {

    @Test
    public void canInterpretDate() {
        assertEquals(WidePatientReader.createInterpretDate("18-apr-2019").toString(), "2019-04-18");
        assertEquals(WidePatientReader.createInterpretDate("17-okt-2018").toString(), "2018-10-17");
    }

    @Test
    public void canLoadEmptyPatient() {

        WideEcrfModel wideEcrfModel;

        wideEcrfModel = ImmutableWideEcrfModel.builder()
                .preTreatments(Lists.newArrayList())
                .biopsies(Lists.newArrayList())
                .treatments(Lists.newArrayList())
                .responses(Lists.newArrayList())
                .build();

        WidePatientReader patientReader = new WidePatientReader(wideEcrfModel, TestCuratorFactory.tumorLocationCurator(),
                TestCuratorFactory.biopsySiteCurator(),
                TestCuratorFactory.treatmentCurator());

        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read("ID", "Melanoma", Lists.newArrayList(sample));

        assertNotNull(patient);
    }
}