package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HospitalModelTest {

    @Test
    public void canReadPIsForCPCTAndDrup() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        final HospitalData hospital = hospitalModel.hospitalPerId("01");
        assertNotNull(hospital);
        assertEquals("CpctPI", HospitalModel.getPI("CPCT02010001", hospital, "henk"));
        assertEquals("DrupPI", HospitalModel.getPI("DRUP01010001", hospital, "henk"));

        // HospitalModelFactory with '*' for drup pi & recipients
        final HospitalData hospital2 = hospitalModel.hospitalPerId("02");
        assertNotNull(hospital2);
        assertEquals("CpctPI2", HospitalModel.getPI("CPCT02010001", hospital2, "henk"));
        assertEquals("CpctPI2", HospitalModel.getPI("DRUP01010001", hospital2, "henk"));

        assertNull(hospitalModel.hospitalPerId("03"));
    }

    @Test
    public void canReadAddress() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        final HospitalData hospital = hospitalModel.hospitalPerId("01");
        assertNotNull(hospital);
        assertEquals("Address", hospital.addressName());
        assertEquals("Zip", hospital.addressZip());
        assertEquals("City", hospital.addressCity());
    }

    @Test
    public void canLookupAddresseeForSample() {
        final HospitalModel hospitalModel = buildTestHospitalModel();

        assertEquals("CpctPI, Address, Zip City", hospitalModel.addresseeStringForSample("henk", "CPCT02010001T"));
    }

    @Test
    public void canLookAddressForCORESample() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("Address, Zip City", hospitalModel.addresseeStringForSample("henk", "CORE18001224T"));
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel() {
        Map<String, HospitalData> hospitalPerId = Maps.newHashMap();
        Map<String, HospitalData> hospitalPerHospitalName = Maps.newHashMap();
        Map<String, HospitalSampleMapping> hospitalPerIdManual = Maps.newHashMap();

        hospitalPerId.put("01", ImmutableHospitalData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));
        hospitalPerId.put("02", ImmutableHospitalData.of("CPCT Recip2", "*", "Address2", "Zip2", "City2", "CpctPI2", "*"));
        hospitalPerHospitalName.put("HMF", ImmutableHospitalData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));
        hospitalPerIdManual.put("CORE18001224T", ImmutableHospitalSampleMapping.of("HMF"));

        return ImmutableHospitalModel.of(hospitalPerId, hospitalPerHospitalName, hospitalPerIdManual);
    }
}
