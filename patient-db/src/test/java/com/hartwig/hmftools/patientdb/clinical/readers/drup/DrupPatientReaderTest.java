package com.hartwig.hmftools.patientdb.clinical.readers.drup;

import static com.hartwig.hmftools.patientdb.clinical.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertNotNull;

import java.time.LocalDate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.clinical.data.Patient;
import com.hartwig.hmftools.patientdb.clinical.data.SampleData;

import org.junit.Test;

public class DrupPatientReaderTest {

    @Test
    public void canReadEmptyPatient() {
        DrupPatientReader patientReader = new DrupPatientReader(
                TestCuratorFactory.primaryTumorCurator(),
                TestCuratorFactory.biopsySiteCurator());

        EcrfPatient ecrfPatient = new EcrfPatient("empty", Maps.newHashMap(), Lists.newArrayList());
        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read(ecrfPatient, Lists.newArrayList(sample));

        assertNotNull(patient);
    }
}