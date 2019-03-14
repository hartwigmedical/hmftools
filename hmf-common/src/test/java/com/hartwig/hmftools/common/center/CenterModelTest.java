package com.hartwig.hmftools.common.center;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CenterModelTest {

    @Test
    public void canReadPIsForCPCTAndDrup() {
        final CenterModel centerModel = buildTestCenterModel();
        final CenterData center = centerModel.centerPerId("01");
        assertNotNull(center);
        assertEquals("CpctPI", CenterModel.getPI("CPCT02010001", center, "henk"));
        assertEquals("DrupPI", CenterModel.getPI("DRUP01010001", center, "henk"));

        // CenterModelFactory with '*' for drup pi & recipients
        final CenterData center2 = centerModel.centerPerId("02");
        assertNotNull(center2);
        assertEquals("CpctPI2", CenterModel.getPI("CPCT02010001", center2, "henk"));
        assertEquals("CpctPI2", CenterModel.getPI("DRUP01010001", center2, "henk"));

        assertNull(centerModel.centerPerId("03"));
    }

    @Test
    public void canReadAddress() {
        final CenterModel centerModel = buildTestCenterModel();
        final CenterData center = centerModel.centerPerId("01");
        assertNotNull(center);
        assertEquals("Address", center.addressName());
        assertEquals("Zip", center.addressZip());
        assertEquals("City", center.addressCity());
    }

    @Test
    public void canLookupAddresseeForSample() {
        final CenterModel centerModel = buildTestCenterModel();

        assertEquals("CpctPI, Address, Zip City", centerModel.addresseeStringForSample("henk", "CPCT02010001T"));
    }

    @Test
    public void canLookAdressForCORESample() {
        final CenterModel centerModel = buildTestCenterModel();
        assertEquals("Address, Zip City", centerModel.addresseeStringForSample("henk", "CORE18001224T"));
    }

    @NotNull
    private static CenterModel buildTestCenterModel() {
        Map<String, CenterData> centerPerId = Maps.newHashMap();
        Map<String, CenterData> centerPerHospital = Maps.newHashMap();
        Map<String, CenterDataManualMapping> centerPerIdManual = Maps.newHashMap();

        centerPerId.put("01", ImmutableCenterData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));
        centerPerId.put("02", ImmutableCenterData.of("CPCT Recip2", "*", "Address2", "Zip2", "City2", "CpctPI2", "*"));
        centerPerHospital.put("HMF", ImmutableCenterData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));
        centerPerIdManual.put("CORE18001224T", ImmutableCenterDataManualMapping.of("HMF"));

        return ImmutableCenterModel.of(centerPerId, centerPerHospital, centerPerIdManual);
    }
}
