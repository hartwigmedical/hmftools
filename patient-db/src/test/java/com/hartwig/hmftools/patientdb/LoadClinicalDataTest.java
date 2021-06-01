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
import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.clinical.validators.PatientValidator;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class LoadClinicalDataTest {

    private static final String TEST_ECRF = Resources.getResource("ecrf/cpct_ecrf.xml").getPath();

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @Test
    public void canLoadUpRealCpctEcrf() throws IOException, XMLStreamException {
        PrimaryTumorCurator primaryTumorCurator = TestCuratorFactory.primaryTumorCurator();
        BiopsySiteCurator biopsySiteCurator = TestCuratorFactory.biopsySiteCurator();
        TreatmentCurator treatmentCurator = TestCuratorFactory.treatmentCurator();

        EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLNoFormStates(TEST_ECRF);
        assertEquals(1, cpctEcrfModel.patientCount());
        assertEquals(1100, Lists.newArrayList(cpctEcrfModel.fields()).size());

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
        assertEquals("CPCT02252500", patient.patientIdentifier());
        assertEquals(0, patient.sequencedBiopsies().size());

        BaselineData baselineData = patient.baselineData();
        assertNotNull(baselineData);
        assertEquals(new Integer(1963), baselineData.birthYear());
        assertEquals("Breast cancer", baselineData.curatedPrimaryTumor().searchTerm());
        assertEquals("female", baselineData.gender());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), baselineData.informedConsentDate());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), baselineData.registrationDate());
        assertEquals(LocalDate.parse("2012-06-22", DATE_FORMATTER), baselineData.deathDate());
        assertEquals("Bernhoven uden", baselineData.hospital());

        PreTreatmentData preTreatmentData = patient.preTreatmentData();
        assertNotNull(preTreatmentData);
        assertEquals("Yes", preTreatmentData.radiotherapyGiven());
        assertEquals("Yes", preTreatmentData.treatmentGiven());
        assertEquals(6, preTreatmentData.drugs().size());
        assertEquals("Bevacizumab", preTreatmentData.drugs().get(0).name());

        List<TumorMarkerData> tumorMarkers = patient.tumorMarkers();
        assertEquals(0, tumorMarkers.size());

        List<BiopsyData> biopsies = patient.clinicalBiopsies();
        assertEquals(1, biopsies.size());
        assertNull(biopsies.get(0).biopsyEvaluable());
        assertEquals("Yes", biopsies.get(0).biopsyTaken());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), biopsies.get(0).date());
        assertEquals("Soft tissue", biopsies.get(0).site());
        assertEquals("near right scapula", biopsies.get(0).location());

        List<BiopsyTreatmentData> treatments = patient.treatments();
        assertEquals(1, treatments.size());
        assertEquals("Yes", treatments.get(0).treatmentGiven());
        assertEquals("No", treatments.get(0).radiotherapyGiven());
        assertEquals(0, treatments.get(0).drugs().size());

        List<BiopsyTreatmentResponseData> responses = patient.treatmentResponses();
        assertEquals(4, responses.size());
        assertNull(responses.get(0).measurementDone());
        assertNull(responses.get(0).boneOnlyDisease());
        assertNull(responses.get(0).responseDate());
        assertEquals(LocalDate.parse("2012-03-05", DATE_FORMATTER), responses.get(1).assessmentDate());
        assertNull(responses.get(0).response());
    }
}
