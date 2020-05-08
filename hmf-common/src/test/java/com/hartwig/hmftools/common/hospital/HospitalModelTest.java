package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HospitalModelTest {

    @Test
    public void extractPIName() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("Someone", hospitalModel.extractHospitalPI("CPCT01010001"));
        assertEquals("Someone", hospitalModel.extractHospitalPI("DRUP01010001"));
        assertEquals("Someone", hospitalModel.extractHospitalPI("WIDE01010001"));
    }

    @Test
    public void extractRequestName() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("Someone1", hospitalModel.extractRequestName("WIDE01010001", "BB"));
        assertNull(hospitalModel.extractRequestName("DRUP01010001", "BB"));
        assertNull(hospitalModel.extractRequestName("CPCT01010001", "BB"));
        assertEquals("BB", hospitalModel.extractRequestName("CORE01010001", "BB"));
    }

    @Test
    public void extractRequestEmail() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("my@email.com", hospitalModel.extractRequestEmail("WIDE01010001", "AA"));
        assertNull(hospitalModel.extractRequestEmail("CPCT01010001", "AA"));
        assertNull(hospitalModel.extractRequestEmail("DRUP01010001", "AA"));
        assertEquals("AA", hospitalModel.extractRequestEmail("CORE01010001", "AA"));
    }

    @Test
    public void canLookupAddresseeForSample() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("Someone, Ext-HMF, 1000 AB AMSTERDAM", hospitalModel.extractHospitalAdress("CPCT02010001T"));
        assertEquals("Someone, Ext-HMF, 1000 AB AMSTERDAM", hospitalModel.extractHospitalAdress("DRUP02010001T"));
        assertEquals("Someone, Ext-HMF, 1000 AB AMSTERDAM", hospitalModel.extractHospitalAdress("WIDE02010001T"));
        assertEquals("Ext-HMF, 1000 AB AMSTERDAM", hospitalModel.extractHospitalAdress("CORE02010001T"));

    }

    @Test
    public void canLookupHospitalNameForSample() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        assertEquals("Ext-HMF", hospitalModel.extractHospitalName("CPCT02010001T"));
        assertEquals("Ext-HMF", hospitalModel.extractHospitalName("DRUP02010001T"));
        assertEquals("Ext-HMF", hospitalModel.extractHospitalName("WIDE02010001T"));
        assertEquals("Ext-HMF", hospitalModel.extractHospitalName("CORE02010001T"));
    }

    @NotNull
    public static HospitalModel buildTestHospitalModel() {

        Map<String, HospitalSampleMapping> sampleHospitalMapping = Maps.newHashMap();
        Map<String, HospitalData> hospitalDataCPCT = Maps.newHashMap();
        Map<String, HospitalData> hospitalDataDRUP = Maps.newHashMap();
        Map<String, HospitalData> hospitalDataWIDE = Maps.newHashMap();
        Map<String, HospitalAdress> hospitalAdress = Maps.newHashMap();

        sampleHospitalMapping.put("CORE18001224T", ImmutableHospitalSampleMapping.of("HOSP1"));
        hospitalDataCPCT.put("01", ImmutableHospitalData.of("01", "Someone", "", ""));
        hospitalDataDRUP.put("01", ImmutableHospitalData.of("01", "Someone", "", ""));
        hospitalDataWIDE.put("01", ImmutableHospitalData.of("01", "Someone", "Someone1", "my@email.com"));
        hospitalAdress.put("01", ImmutableHospitalAdress.of("01", "Ext-HMF","1000 AB","AMSTERDAM"));

        return ImmutableHospitalModel.of(
                sampleHospitalMapping,
                hospitalDataCPCT,
                hospitalDataDRUP,
                hospitalDataWIDE,
                hospitalAdress);
    }
}
