package com.hartwig.hmftools.patientdb.readers.drup;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertNotNull;

import java.time.LocalDate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.junit.Test;

public class DrupPatientReaderTest {

    @Test
    public void canReadEmptyPatient() {
        DrupPatientReader patientReader =
                new DrupPatientReader(TestCuratorFactory.tumorLocationCurator(), TestCuratorFactory.biopsySiteCurator());

        EcrfPatient ecrfPatient = new EcrfPatient("empty", Maps.newHashMap(), Lists.newArrayList());
        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read(ecrfPatient, Lists.newArrayList(sample));

        assertNotNull(patient);
    }
}