package com.hartwig.hmftools.patientdb.clinical.readers.drup;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertNotNull;

import java.time.LocalDate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientdb.clinical.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;

import org.junit.Test;

public class DrupPatientReaderTest {

    private static final String INFORMED_CONSENTS_TSV = Resources.getResource("consents/informed_consents.tsv").getPath();

    @Test
    public void canReadEmptyPatient() {
        DrupPatientReader patientReader = new DrupPatientReader(
                TestCuratorFactory.primaryTumorCurator(),
                TestCuratorFactory.biopsySiteCurator());

        EcrfPatient ecrfPatient = new EcrfPatient("empty", Maps.newHashMap(), Lists.newArrayList());
        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read(ecrfPatient, Lists.newArrayList(sample), INFORMED_CONSENTS_TSV);

        assertNotNull(patient);
    }
}