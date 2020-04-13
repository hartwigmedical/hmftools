package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.*;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWidePreTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideResponseData;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WidePatientReader;
import com.hartwig.hmftools.patientdb.readers.wide.WidePreTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;
import com.hartwig.hmftools.patientdb.readers.wide.WideTreatmentData;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class WidePatientReaderTest {

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

    @Test
    public void canReadCombiPreTreatment() {
        List<WidePreTreatmentData> preTreatmentDataCombi = Lists.newArrayList(ImmutableWidePreTreatmentData.builder()
                .patientId("WIDE00000001")
                .previousTherapy("1")
                .drug1("carboplatin")
                .drug2("erlotinib")
                .drug3(Strings.EMPTY)
                .drug4(Strings.EMPTY)
                .dateLastSystemicTherapy("20-nov-2015")
                .build());
        List<WidePreTreatmentData> preTreatmentData = Lists.newArrayList(ImmutableWidePreTreatmentData.builder()
                .patientId("WIDE00000001")
                .previousTherapy("1")
                .drug1("carboplatin")
                .drug2(Strings.EMPTY)
                .drug3(Strings.EMPTY)
                .drug4(Strings.EMPTY)
                .dateLastSystemicTherapy("20-nov-2015")
                .build());
        TreatmentCurator treatmentCurator = TestCuratorFactory.treatmentCurator();
        String patientIdentifier = "WIDE00000001";
        LocalDate biopsyDate = WidePatientReader.createInterpretDateEN("21-May-2019");
        List<WideTreatmentData> treatmentData = Lists.newArrayList();

        assertEquals(WidePatientReader.readDrugsPreTreatment(preTreatmentDataCombi,
                treatmentCurator,
                patientIdentifier,
                biopsyDate,
                treatmentData).get(0).name(), "carboplatin,erlotinib");

        assertEquals(WidePatientReader.readDrugsPreTreatment(preTreatmentData,
                treatmentCurator,
                patientIdentifier,
                biopsyDate,
                treatmentData).get(0).name(), "carboplatin");

    }

    @Test
    public void canLoadEmptyPatient() {
        Lims lims = LimsFactory.empty();
        WideEcrfModel wideEcrfModel;

        wideEcrfModel = ImmutableWideEcrfModel.builder()
                .preTreatments(Lists.newArrayList())
                .biopsies(Lists.newArrayList())
                .treatments(Lists.newArrayList())
                .responses(Lists.newArrayList())
                .build();

        WidePatientReader patientReader = new WidePatientReader(wideEcrfModel,
                TestCuratorFactory.tumorLocationCurator(),
                TestCuratorFactory.biopsySiteCurator(),
                TestCuratorFactory.treatmentCurator());

        SampleData sample = sampleBuilder(LocalDate.parse("2017-01-01")).build();

        Patient patient = patientReader.read("ID", "Melanoma", Lists.newArrayList(sample), lims, "FR123");

        assertNotNull(patient);
    }
}