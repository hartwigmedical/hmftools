package com.hartwig.hmftools.common.center;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CenterModelTest {

    @Test
    public void canReadPIsForCPCTAndDrup() {
        final CenterModel centerModel = buildTestCenterModel();
        final CenterData center = centerModel.centerPerId("01");
        assertNotNull(center);
        assertEquals("CpctPI", CenterModel.getPI("CPCT02010001", center));
        assertEquals("DrupPI", CenterModel.getPI("DRUP01010001", center));
        assertEquals(Strings.EMPTY, CenterModel.getPI("DoesNotExist", center));

        // Center with '*' for drup pi & recipients
        final CenterData center2 = centerModel.centerPerId("02");
        assertNotNull(center2);
        assertEquals("CpctPI2", CenterModel.getPI("CPCT02010001", center2));
        assertEquals("CpctPI2", CenterModel.getPI("DRUP01010001", center2));
        assertEquals(Strings.EMPTY, CenterModel.getPI("DoesNotExist", center2));

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

        assertEquals("CpctPI, Address, Zip City", centerModel.addresseeStringForSample("CPCT02010001T"));
        assertNull(centerModel.addresseeStringForSample("DoesNotExist"));
    }

    @Test
    public void canLookupAddresseeForProject() {
        final CenterModel centerModel = buildTestCenterModel();

        assertEquals("Address, Zip City", centerModel.addresseeStringForProject("HMF-001-002"));
        assertNull(centerModel.addresseeStringForProject("DoesNotExist"));
    }

    @NotNull
    private static CenterModel buildTestCenterModel() {
        Map<String, CenterData> centerPerId = Maps.newHashMap();
        Map<String, CenterData> centerPerHospital = Maps.newHashMap();

        centerPerId.put("01", ImmutableCenterData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));
        centerPerId.put("02", ImmutableCenterData.of("CPCT Recip2", "*", "Address2", "Zip2", "City2", "CpctPI2", "*"));
        centerPerHospital.put("HMF", ImmutableCenterData.of("CPCT Recip", "DRUP Recip", "Address", "Zip", "City", "CpctPI", "DrupPI"));

        return ImmutableCenterModel.of(centerPerId, centerPerHospital);
    }
}
