package com.hartwig.hmftools.patientdb.readers.cpct;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PatientReaderTest {

    private static final String TEST_ECRF = Resources.getResource("test_ecrf.xml").getPath();

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @Test
    public void canReadCpctBaselineInfo() {
        final EcrfModel model = loadTestEcrf();
        assertEquals(1, model.patientCount());

        final Map<Integer, String> hospitals = Util.extractHospitalMap(model);
        final BaselineReader baselineReader = new BaselineReader(TestCuratorFactory.tumorLocationCurator(), hospitals);

        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final BaselineData baselineData = baselineReader.read(cpctPatient);
        assertEquals("Breast cancer", baselineData.curatedTumorLocation().searchTerm());
        assertEquals("Breast", baselineData.curatedTumorLocation().primaryTumorLocation());
        assertEquals("Subtype unknown", baselineData.curatedTumorLocation().subType());
        assertEquals("female", baselineData.gender());
        assertEquals("Bernhoven uden", baselineData.hospital());
        assertEquals(new Integer(1963), baselineData.birthYear());
        assertEquals(LocalDate.parse("2012-06-22", DATE_FORMATTER), baselineData.deathDate());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), baselineData.registrationDate());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), baselineData.informedConsentDate());
    }

    @Test
    public void canReadCpctPatientPreTherapy() throws IOException {
        final EcrfModel model = loadTestEcrf();
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final PreTreatmentData preTreatmentData = new PreTreatmentReader(TestCuratorFactory.treatmentCurator()).read(cpctPatient);

        assertEquals("Yes", preTreatmentData.treatmentGiven());
        assertEquals("Yes", preTreatmentData.radiotherapyGiven());
        final List<DrugData> drugs = preTreatmentData.drugs();
        assertEquals(6, drugs.size());
    }

    @Test
    public void canReadCpctPatientTreatments() throws IOException {
        final EcrfModel model = loadTestEcrf();
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final List<BiopsyTreatmentData> treatments = new BiopsyTreatmentReader(TestCuratorFactory.treatmentCurator()).read(cpctPatient);

        assertEquals(1, treatments.size());
        assertEquals(1, treatments.get(0).drugs().size());
        final LocalDate startDate = LocalDate.parse("2012-02-18", DATE_FORMATTER);
        final LocalDate endDate = LocalDate.parse("2012-04-02", DATE_FORMATTER);
        assertEquals("Bevacizumab", treatments.get(0).drugs().get(0).name());
        assertEquals(startDate, treatments.get(0).drugs().get(0).startDate());
        assertEquals(endDate, treatments.get(0).drugs().get(0).endDate());
        assertEquals(startDate, treatments.get(0).startDate());
        assertEquals(endDate, treatments.get(0).endDate());
        assertEquals("Bevacizumab", treatments.get(0).treatmentName());
        assertEquals("Yes", treatments.get(0).treatmentGiven());
    }

    @Test
    public void canReadCpctPatientTreatmentResponses() {
        final EcrfModel model = loadTestEcrf();
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(cpctPatient);

        assertEquals(4, treatmentResponses.size());
        assertEquals(LocalDate.parse("2012-02-09", DATE_FORMATTER), treatmentResponses.get(0).date());
        assertEquals(null, treatmentResponses.get(0).response());
        assertEquals(null, treatmentResponses.get(0).measurementDone());

        assertEquals(LocalDate.parse("2012-03-05", DATE_FORMATTER), treatmentResponses.get(1).date());
        assertEquals("PD", treatmentResponses.get(1).response());
        assertEquals(null, treatmentResponses.get(1).measurementDone());

        assertEquals(LocalDate.parse("2012-04-23", DATE_FORMATTER), treatmentResponses.get(2).date());
        assertEquals("SD", treatmentResponses.get(2).response());
        assertEquals(null, treatmentResponses.get(2).measurementDone());

        assertEquals(LocalDate.parse("2012-06-08", DATE_FORMATTER), treatmentResponses.get(3).date());
        assertEquals("PR", treatmentResponses.get(3).response());
        assertEquals("Yes", treatmentResponses.get(3).measurementDone());
    }

    @Test
    public void canReadTumorMarkers() {
        final EcrfModel model = loadTestEcrf();
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final List<TumorMarkerData> tumorMarkers = TumorMarkerReader.read(cpctPatient);

        assertEquals(0, tumorMarkers.size());
    }

    @NotNull
    private static EcrfModel loadTestEcrf() {
        EcrfModel model = null;
        try {
            model = EcrfModel.loadFromXMLWithFormStates(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
        } catch (XMLStreamException | FileNotFoundException e) {
            // KODU: Ignore, should not happen in test.
        }

        if (model == null) {
            throw new IllegalStateException("No cpct ecrf model could be constructed!");
        }

        return model;
    }
}
