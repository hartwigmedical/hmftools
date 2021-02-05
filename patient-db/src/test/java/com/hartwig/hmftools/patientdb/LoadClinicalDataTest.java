package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.TumorMarkerData;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.clinical.validators.PatientValidator;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class LoadClinicalDataTest {

    private static final String TEST_ECRF = Resources.getResource("test_cpct_ecrf.xml").getPath();

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @Test
    public void canLoadUpRealCpctEcrf() throws IOException, XMLStreamException {
        PrimaryTumorCurator primaryTumorCurator = TestCuratorFactory.primaryTumorCurator();
        BiopsySiteCurator biopsySiteCurator = TestCuratorFactory.biopsySiteCurator();
        TreatmentCurator treatmentCurator = TestCuratorFactory.treatmentCurator();

        EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLNoFormStates(TEST_ECRF);
        assertEquals(1, cpctEcrfModel.patientCount());
        assertEquals(1298, Lists.newArrayList(cpctEcrfModel.fields()).size());

        EcrfPatientReader cpctPatientReader = new CpctPatientReader(primaryTumorCurator,
                CpctUtil.extractHospitalMap(cpctEcrfModel),
                biopsySiteCurator,
                treatmentCurator);

        List<ValidationFinding> allFindings = Lists.newArrayList();
        for (EcrfPatient ecrfPatient : cpctEcrfModel.patients()) {
            Patient patient = cpctPatientReader.read(ecrfPatient, Lists.newArrayList());
            assertPatient(patient);
            allFindings.addAll(PatientValidator.validatePatient(patient));
        }

        assertTrue(allFindings.size() > 0);
    }

    private static void assertPatient(@Nullable Patient patient) {
        assertNotNull(patient);
        assertEquals("CPCT02000015", patient.patientIdentifier());
        assertEquals(0, patient.sequencedBiopsies().size());

        BaselineData baselineData = patient.baselineData();
        assertNotNull(baselineData);
        assertEquals(new Integer(1984), baselineData.birthYear());
        assertEquals("Gastrointestinal Stromal Tumors (GIST)", baselineData.curatedPrimaryTumor().searchTerm());
        assertEquals("female", baselineData.gender());
        assertEquals(LocalDate.parse("2018-06-02", DATE_FORMATTER), baselineData.informedConsentDate());
        assertEquals(LocalDate.parse("2018-12-10", DATE_FORMATTER), baselineData.registrationDate());
        assertNull(baselineData.deathDate());
        assertNull(baselineData.hospital());

        PreTreatmentData preTreatmentData = patient.preTreatmentData();
        assertNotNull(preTreatmentData);
        assertEquals("No", preTreatmentData.radiotherapyGiven());
        assertEquals("Yes", preTreatmentData.treatmentGiven());
        assertEquals(1, preTreatmentData.drugs().size());
        assertEquals("Imatinib", preTreatmentData.drugs().get(0).name());

        List<TumorMarkerData> tumorMarkers = patient.tumorMarkers();
        assertEquals(0, tumorMarkers.size());

        List<BiopsyData> biopsies = patient.clinicalBiopsies();
        assertEquals(1, biopsies.size());
        assertNull(biopsies.get(0).biopsyEvaluable());
        assertEquals("Yes", biopsies.get(0).biopsyTaken());
        assertEquals(LocalDate.parse("2017-07-14", DATE_FORMATTER), biopsies.get(0).date());
        assertEquals("abdomen", biopsies.get(0).site());
        assertEquals("left mesenterial", biopsies.get(0).location());

        List<BiopsyTreatmentData> treatments = patient.treatments();
        assertEquals(1, treatments.size());
        assertEquals("Yes", treatments.get(0).treatmentGiven());
        assertEquals("No", treatments.get(0).radiotherapyGiven());
        assertEquals(1, treatments.get(0).drugs().size());
        assertEquals(LocalDate.parse("2014-11-16", DATE_FORMATTER), treatments.get(0).drugs().get(0).startDate());
        assertNull(treatments.get(0).drugs().get(0).endDate());
        assertEquals("Imatinib", treatments.get(0).drugs().get(0).name());

        List<BiopsyTreatmentResponseData> responses = patient.treatmentResponses();
        assertEquals(7, responses.size());
        assertEquals("Yes", responses.get(0).measurementDone());
        assertNull(responses.get(0).boneOnlyDisease());
        assertEquals(LocalDate.parse("2012-05-21", DATE_FORMATTER), responses.get(0).responseDate());
        assertEquals(LocalDate.parse("2015-08-30", DATE_FORMATTER), responses.get(1).assessmentDate());
        assertEquals("SD", responses.get(0).response());
    }
}
