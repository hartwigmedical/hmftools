package com.hartwig.hmftools.patientdb.clinical.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.lims.Lims;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HospitalModelTest {

    private static final String HOSPITAL_NAME = "HMF Name";
    private static final String HOSPITAL_ZIP = "HMF Zip";
    private static final String HOSPITAL_CITY = "HMF City";
    private static final String EXPECTED_HOSPITAL_ADDRESS = HOSPITAL_ZIP + " " + HOSPITAL_CITY;

    private static final String CPCT_PI = "CPCT-PI";
    private static final String DRUP_PI = "DRUP-PI";
    private static final String WIDE_PI = "WIDE-PI";

    private static final String WIDE_REQUESTER_NAME = "WIDE Requester";
    private static final String WIDE_REQUESTER_EMAIL = "wide@wide.com";
    private static final String CORE_REQUESTER_NAME = "CORE Requester";
    private static final String CORE_REQUESTER_EMAIL = "core@core.com";

    private static final String MANUALLY_MAPPED_CORE_SAMPLE = "COREManualMap";

    @Test
    public void canExtractHospitalDataFromHospitalModel() {
        String hospitalId = "01";
        HospitalModel model = buildTestHospitalModel(hospitalId);

        LimsCohortConfig cohortConfigCPCT = LimsCohortTestFactory.createConfigForHospitalModel("CPCT", false, false);
        HospitalContactData dataCPCT =
                model.queryHospitalData("CPCT02010001T", cohortConfigCPCT, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertEquals(CPCT_PI, dataCPCT.hospitalPI());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCPCT.requesterName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCPCT.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCPCT.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCPCT.hospitalAddress());

        LimsCohortConfig cohortConfigDRUP = LimsCohortTestFactory.createConfigForHospitalModel("DRUP", false, false);
        HospitalContactData dataDRUP =
                model.queryHospitalData("DRUP02010001T", cohortConfigDRUP, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertEquals(DRUP_PI, dataDRUP.hospitalPI());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataDRUP.requesterName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataDRUP.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataDRUP.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataDRUP.hospitalAddress());

        LimsCohortConfig cohortConfigWIDE = LimsCohortTestFactory.createConfigForHospitalModel("WIDE", false, false);
        HospitalContactData dataWIDE =
                model.queryHospitalData("WIDE02010001T", cohortConfigWIDE, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertEquals(WIDE_PI, dataWIDE.hospitalPI());
        assertEquals(WIDE_REQUESTER_NAME, dataWIDE.requesterName());
        assertEquals(WIDE_REQUESTER_EMAIL, dataWIDE.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataWIDE.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataWIDE.hospitalAddress());

        LimsCohortConfig cohortConfigCORE = LimsCohortTestFactory.createConfigForHospitalModel("CORE", false, true);
        HospitalContactData dataCORE =
                model.queryHospitalData("CORE02010001T", cohortConfigCORE, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCORE.hospitalPI());
        assertEquals(CORE_REQUESTER_NAME, dataCORE.requesterName());
        assertEquals(CORE_REQUESTER_EMAIL, dataCORE.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCORE.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCORE.hospitalAddress());

        HospitalContactData dataCOREManuallyMapped =
                model.queryHospitalData(MANUALLY_MAPPED_CORE_SAMPLE, cohortConfigCORE, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCOREManuallyMapped.hospitalPI());
        assertEquals(CORE_REQUESTER_NAME, dataCOREManuallyMapped.requesterName());
        assertEquals(CORE_REQUESTER_EMAIL, dataCOREManuallyMapped.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCOREManuallyMapped.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCOREManuallyMapped.hospitalAddress());

        HospitalContactData dataSampleDoesNotExist =
                model.queryHospitalData("I Don't exist", cohortConfigCORE, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertNull(dataSampleDoesNotExist);

        HospitalContactData dataHospitalDoesNotExist =
                model.queryHospitalData("CPCT02020202T", cohortConfigCPCT, CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL);
        assertNull(dataHospitalDoesNotExist);
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel(@NotNull String hospitalId) {
        Map<String, HospitalAddress> hospitalAddress = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactCPCT = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactDRUP = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactWIDE = Maps.newHashMap();
        Map<String, String> sampleHospitalMapping = Maps.newHashMap();

        hospitalAddress.put(hospitalId,
                ImmutableHospitalAddress.builder()
                        .hospitalName(HOSPITAL_NAME)
                        .hospitalZip(HOSPITAL_ZIP)
                        .hospitalCity(HOSPITAL_CITY)
                        .build());
        hospitalContactCPCT.put(hospitalId, ImmutableHospitalPersons.builder().hospitalPI(CPCT_PI).build());
        hospitalContactDRUP.put(hospitalId, ImmutableHospitalPersons.builder().hospitalPI(DRUP_PI).build());
        hospitalContactWIDE.put(hospitalId,
                ImmutableHospitalPersons.builder()
                        .hospitalPI(WIDE_PI)
                        .requesterName(WIDE_REQUESTER_NAME)
                        .requesterEmail(WIDE_REQUESTER_EMAIL)
                        .build());
        sampleHospitalMapping.put(MANUALLY_MAPPED_CORE_SAMPLE, hospitalId);

        return ImmutableHospitalModel.builder()
                .hospitalAddressMap(hospitalAddress)
                .hospitalPersonsCPCT(hospitalContactCPCT)
                .hospitalPersonsDRUP(hospitalContactDRUP)
                .hospitalPersonsWIDE(hospitalContactWIDE)
                .sampleToHospitalMapping(sampleHospitalMapping)
                .build();
    }
}
