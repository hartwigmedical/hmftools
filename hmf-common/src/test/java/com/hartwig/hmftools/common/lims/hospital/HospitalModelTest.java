package com.hartwig.hmftools.common.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
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
    @Ignore
    public void canExtractHospitalDataFromHospitalModel() {
        String hospitalId = "01";
        HospitalModel model = buildTestHospitalModel(hospitalId);

        HospitalContactData dataCPCT =
                model.queryHospitalData("CPCT02010001T", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, buildTestCohortModel("CPCT").queryCohortData("CPCT"));
        assertEquals(CPCT_PI, dataCPCT.hospitalPI());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCPCT.requesterName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCPCT.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCPCT.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCPCT.hospitalAddress());

        HospitalContactData dataDRUP =
                model.queryHospitalData("DRUP02010001T", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, buildTestCohortModel("DRUP").queryCohortData("DRUP"));
        assertEquals(DRUP_PI, dataDRUP.hospitalPI());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataDRUP.requesterName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataDRUP.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataDRUP.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataDRUP.hospitalAddress());

        HospitalContactData dataWIDE =
                model.queryHospitalData("WIDE02010001T", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, buildTestCohortModel("WIDE").queryCohortData("WIDE"));
        assertEquals(WIDE_PI, dataWIDE.hospitalPI());
        assertEquals(WIDE_REQUESTER_NAME, dataWIDE.requesterName());
        assertEquals(WIDE_REQUESTER_EMAIL, dataWIDE.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataWIDE.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataWIDE.hospitalAddress());

        HospitalContactData dataCORE =
                model.queryHospitalData("CORE02010001T", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, buildTestCohortModel("CORE").queryCohortData("CORE"));
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCORE.hospitalPI());
        assertEquals(CORE_REQUESTER_NAME, dataCORE.requesterName());
        assertEquals(CORE_REQUESTER_EMAIL, dataCORE.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCORE.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCORE.hospitalAddress());

        HospitalContactData dataCOREManuallyMapped = model.queryHospitalData(MANUALLY_MAPPED_CORE_SAMPLE,
                CORE_REQUESTER_NAME,
                CORE_REQUESTER_EMAIL,
                buildTestCohortModel("CORE").queryCohortData("CORE"));
        assertEquals(Lims.NOT_AVAILABLE_STRING, dataCOREManuallyMapped.hospitalPI());
        assertEquals(CORE_REQUESTER_NAME, dataCOREManuallyMapped.requesterName());
        assertEquals(CORE_REQUESTER_EMAIL, dataCOREManuallyMapped.requesterEmail());
        assertEquals(HOSPITAL_NAME, dataCOREManuallyMapped.hospitalName());
        assertEquals(EXPECTED_HOSPITAL_ADDRESS, dataCOREManuallyMapped.hospitalAddress());

        HospitalContactData dataSampleDoesNotExist =
                model.queryHospitalData("I Don't exist", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, null);
        assertNull(dataSampleDoesNotExist);

        HospitalContactData dataHospitalDoesNotExist =
                model.queryHospitalData("CPCT02020202T", CORE_REQUESTER_NAME, CORE_REQUESTER_EMAIL, buildTestCohortModel("CPCT").queryCohortData("CPCT"));
        assertNull(dataHospitalDoesNotExist);
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel(@NotNull String hospitalId) {
        Map<String, HospitalAddress> hospitalAddress = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactCPCT = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactDRUP = Maps.newHashMap();
        Map<String, HospitalPersons> hospitalContactWIDE = Maps.newHashMap();
        Map<String, String> sampleHospitalMapping = Maps.newHashMap();

        hospitalAddress.put(hospitalId, ImmutableHospitalAddress.of(HOSPITAL_NAME, HOSPITAL_ZIP, HOSPITAL_CITY));
        hospitalContactCPCT.put(hospitalId, ImmutableHospitalPersons.of(CPCT_PI, null, null));
        hospitalContactDRUP.put(hospitalId, ImmutableHospitalPersons.of(DRUP_PI, null, null));
        hospitalContactWIDE.put(hospitalId, ImmutableHospitalPersons.of(WIDE_PI, WIDE_REQUESTER_NAME, WIDE_REQUESTER_EMAIL));
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
