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
    public void canDeterminePIForCPCTAndDrupAndWide() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        final HospitalData hospital = hospitalModel.hospitalPerId("01");
        assertNotNull(hospital);
        assertEquals("CpctPI", HospitalModel.determinePI("CPCT02010001", hospital));
        assertEquals("DrupPI", HospitalModel.determinePI("DRUP01010001", hospital));
        assertEquals("WidePI", HospitalModel.determinePI("WIDE01010001", hospital));

        // Revert to CPCT PI with '*' for DRUP PI & recipients
        final HospitalData hospital2 = hospitalModel.hospitalPerId("02");
        assertNotNull(hospital2);
        assertEquals("CpctPI2", HospitalModel.determinePI("CPCT02010001", hospital2));
        assertEquals("CpctPI2", HospitalModel.determinePI("DRUP01010001", hospital2));

        assertNull(hospitalModel.hospitalPerId("03"));
    }

    @Test
    public void canDetermineEmailPIForCPCTAndDrupAndWide() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        final HospitalData hospital = hospitalModel.hospitalPerId("01");
        assertNotNull(hospital);
        assertEquals("WIDE Recip", HospitalModel.determinePIEmail("WIDE01010001", hospital));
        assertEquals("CPCT Recip", HospitalModel.determinePIEmail("CPCT02010001", hospital));
        assertEquals("DRUP Recip", HospitalModel.determinePIEmail("DRUP01010001", hospital));

        // Revert to CPCT PI with '*' for DRUP PI & recipients
        final HospitalData hospital2 = hospitalModel.hospitalPerId("02");
        assertNotNull(hospital2);
        assertEquals("CPCT Recip2", HospitalModel.determinePIEmail("CPCT02010001", hospital2));
    }

    @Test
    public void extractHospitalName() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("ExtHosp1", hospitalModel.externalHospitalName("WIDE01010001"));
    }

    @Test
    public void extractPIName() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("WidePI", hospitalModel.PIName("WIDE01010001"));
    }

    @Test
    public void extractPIEmail() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("WIDE Recip", hospitalModel.PIEmail("WIDE01010001"));
    }

    @Test
    public void canReadHospitalNameAndAddress() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        final HospitalData hospital = hospitalModel.hospitalPerId("01");
        assertNotNull(hospital);
        assertEquals("ExtHosp1", hospital.externalHospitalName());
        assertEquals("Zip", hospital.addressZip());
        assertEquals("City", hospital.addressCity());
    }

    @Test
    public void canLookupAddresseeForSample() {
        final HospitalModel hospitalModel = buildTestHospitalModel();

        assertEquals("CpctPI, ExtHosp1, Zip City", hospitalModel.fullAddresseeString("CPCT02010001T"));
    }

    @Test
    public void canLookupAddressForCORESample() {
        final HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("ExtHosp1, Zip City", hospitalModel.fullAddresseeString("CORE18001224T"));
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel() {
        Map<String, HospitalData> hospitalPerId = Maps.newHashMap();
        Map<String, HospitalSampleMapping> hospitalPerIdManual = Maps.newHashMap();

        hospitalPerId.put("01",
                ImmutableHospitalData.of("HOSP1",
                        "CPCT Recip",
                        "DRUP Recip",
                        "WIDE Recip",
                        "ExtHosp1",
                        "Zip",
                        "City",
                        "CpctPI",
                        "DrupPI",
                        "WidePI"));
        hospitalPerId.put("02",
                ImmutableHospitalData.of("HOSP2", "CPCT Recip2", "*", "", "ExtHosp2", "Zip2", "City2", "CpctPI2", "*", "WidePI2"));
        hospitalPerIdManual.put("CORE18001224T", ImmutableHospitalSampleMapping.of("HOSP1"));

        return ImmutableHospitalModel.of(hospitalPerId, hospitalPerIdManual);
    }
}
