package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.assertEquals;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.PatientData;

import org.junit.Test;

public class PatientReaderTest {
    private static final String TEST_ECRF = Resources.getResource("ecrf.xml").getPath();
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");
    private static final String TREATMENT_MAPPING_CSV = Resources.getResource("treatment_mapping.csv").getPath();
    private static final String TUMOR_LOCATION_MAPPING_CSV = Resources.getResource("tumor_location_mapping.csv").getPath();

    @Test
    public void canReadCpctPatientInfo() throws IOException, XMLStreamException {
        final CpctEcrfModel model = CpctEcrfModel.loadFromXML(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final CpctPatientReader cpctPatientReader =
                new CpctPatientReader(model, new TumorLocationCurator(new FileInputStream(TUMOR_LOCATION_MAPPING_CSV)));
        final PatientData patientData = cpctPatientReader.read(cpctPatient);
        assertEquals("CPCT02252500", patientData.cpctId());
        assertEquals("Breast cancer", patientData.primaryTumorLocation().searchTerm());
        assertEquals("Breast", patientData.primaryTumorLocation().category());
        assertEquals("Breast Cancer: subtype unknown", patientData.primaryTumorLocation().subcategory());
        assertEquals("female", patientData.gender());
        assertEquals("Bernhoven uden", patientData.hospital());
        assertEquals(new Integer(1963), patientData.birthYear());
        assertEquals(LocalDate.parse("2012-06-22", DATE_FORMATTER), patientData.deathDate());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), patientData.registrationDate());
    }

    @Test
    public void canReadCpctPatientBiopsies() throws FileNotFoundException, XMLStreamException {
        final CpctEcrfModel model = CpctEcrfModel.loadFromXML(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final List<BiopsyData> biopsies = BiopsyReader.read(cpctPatient);
        assertEquals(1, biopsies.size());
        assertEquals("Soft tissue", biopsies.get(0).site());
        assertEquals("near right scapula", biopsies.get(0).location());
        assertEquals(LocalDate.parse("2012-02-17", DATE_FORMATTER), biopsies.get(0).date());
    }

    @Test
    public void canReadCpctPatientTreatments() throws IOException, XMLStreamException {
        final CpctEcrfModel model = CpctEcrfModel.loadFromXML(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
        assertEquals(1, model.patientCount());
        final EcrfPatient cpctPatient = model.patients().iterator().next();
        final List<BiopsyTreatmentData> treatments =
                new BiopsyTreatmentReader(new TreatmentCurator(new FileInputStream(TREATMENT_MAPPING_CSV))).read(cpctPatient);
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
    public void canReadCpctPatientTreatmentResponses() throws FileNotFoundException, XMLStreamException {
        final CpctEcrfModel model = CpctEcrfModel.loadFromXML(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
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
}
