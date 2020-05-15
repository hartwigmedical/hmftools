package com.hartwig.hmftools.common.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HospitalModelTest {

    @Test
    public void canExtractHospitalDataFromHospitalModel() {
        HospitalModel model = buildTestHospitalModel();

        HospitalData dataCPCT = model.queryHospitalData("CPCT02010001T", "coreRequester", "coreRequesterEmail");
        assertEquals("CPCT-PI", dataCPCT.hospitalPI());
        assertNull(dataCPCT.requesterName());
        assertNull(dataCPCT.requesterEmail());
        assertEquals("HMF", dataCPCT.hospitalName());
        assertEquals("1000 AB AMSTERDAM", dataCPCT.hospitalAddress());

        HospitalData dataDRUP = model.queryHospitalData("DRUP02010001T", "coreRequester", "coreRequesterEmail");
        assertEquals("DRUP-PI", dataDRUP.hospitalPI());
        assertNull(dataDRUP.requesterName());
        assertNull(dataDRUP.requesterEmail());
        assertEquals("HMF", dataDRUP.hospitalName());
        assertEquals("1000 AB AMSTERDAM", dataDRUP.hospitalAddress());

        HospitalData dataWIDE = model.queryHospitalData("WIDE02010001T", "coreRequester", "coreRequesterEmail");
        assertEquals("WIDE-PI", dataWIDE.hospitalPI());
        assertEquals("WIDE-req", dataWIDE.requesterName());
        assertEquals("wide@email.com", dataWIDE.requesterEmail());
        assertEquals("HMF", dataWIDE.hospitalName());
        assertEquals("1000 AB AMSTERDAM", dataWIDE.hospitalAddress());

        HospitalData dataCORE = model.queryHospitalData("CORE02010001T", "coreRequester", "coreRequesterEmail");
        assertNull(dataCORE.hospitalPI());
        assertEquals("coreRequester", dataCORE.requesterName());
        assertEquals("coreRequesterEmail", dataCORE.requesterEmail());
        assertEquals("HMF", dataCORE.hospitalName());
        assertEquals("1000 AB AMSTERDAM", dataCORE.hospitalAddress());

        HospitalData dataCOREManuallyMapped = model.queryHospitalData("CORE18123456T", "coreRequester", "coreRequesterEmail");
        assertNull(dataCOREManuallyMapped.hospitalPI());
        assertEquals("coreRequester", dataCOREManuallyMapped.requesterName());
        assertEquals("coreRequesterEmail", dataCOREManuallyMapped.requesterEmail());
        assertEquals("HMF", dataCOREManuallyMapped.hospitalName());
        assertEquals("1000 AB AMSTERDAM", dataCOREManuallyMapped.hospitalAddress());

        HospitalData dataSampleDoesNotExist = model.queryHospitalData("I Don't exist", "coreRequester", "coreRequesterEmail");
        assertNull(dataSampleDoesNotExist);

        HospitalData dataHospitalDoesNotExist = model.queryHospitalData("CPCT02020202T", "coreRequester", "coreRequesterEmail");
        assertNull(dataHospitalDoesNotExist);
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel() {
        Map<String, HospitalAddress> hospitalAddress = Maps.newHashMap();
        Map<String, HospitalContact> hospitalContactCPCT = Maps.newHashMap();
        Map<String, HospitalContact> hospitalContactDRUP = Maps.newHashMap();
        Map<String, HospitalContact> hospitalContactWIDE = Maps.newHashMap();
        Map<String, String> sampleHospitalMapping = Maps.newHashMap();

        hospitalAddress.put("01", ImmutableHospitalAddress.of("01", "HMF", "1000 AB", "AMSTERDAM"));
        hospitalContactCPCT.put("01", ImmutableHospitalContact.of("01", "CPCT-PI", null, null));
        hospitalContactDRUP.put("01", ImmutableHospitalContact.of("01", "DRUP-PI", null, null));
        hospitalContactWIDE.put("01", ImmutableHospitalContact.of("01", "WIDE-PI", "WIDE-req", "wide@email.com"));
        sampleHospitalMapping.put("CORE18123456T", "01");

        return ImmutableHospitalModel.builder()
                .hospitalAddress(hospitalAddress)
                .hospitalContactCPCT(hospitalContactCPCT)
                .hospitalContactDRUP(hospitalContactDRUP)
                .hospitalContactWIDE(hospitalContactWIDE)
                .sampleToHospitalMapping(sampleHospitalMapping)
                .build();
    }
}
