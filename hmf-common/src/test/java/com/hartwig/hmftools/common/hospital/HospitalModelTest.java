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
        HospitalModel hospitalModel = buildTestHospitalModel();
       // HospitalData hospital = hospitalModel.hospitalPerId("01");
       // assertNotNull(hospital);
        //        assertEquals("CpctPI", HospitalModel.determinePIName("CPCT02010001", hospital));
        //        assertEquals("DrupPI", HospitalModel.determinePIName("DRUP01010001", hospital));
        //        assertEquals("WidePI", HospitalModel.determinePIName("WIDE01010001", hospital));

        // Revert to CPCT PI with '*' for DRUP PI & recipients
      //  HospitalData hospital2 = hospitalModel.hospitalPerId("02");
       // assertNotNull(hospital2);
        //        assertEquals("CpctPI2", HospitalModel.determinePIName("CPCT02010001", hospital2));
        //        assertEquals("CpctPI2", HospitalModel.determinePIName("DRUP01010001", hospital2));

     //   assertNull(hospitalModel.hospitalPerId("03"));
    }

    @Test
    public void canDetermineEmailPIForCPCTAndDrupAndWide() {
        HospitalModel hospitalModel = buildTestHospitalModel();
      //  HospitalData hospital = hospitalModel.hospitalPerId("01");
      //  assertNotNull(hospital);
        //        assertEquals("WIDE Recip", HospitalModel.determinePIEmail("WIDE01010001", hospital));
        //        assertEquals("CPCT Recip", HospitalModel.determinePIEmail("CPCT02010001", hospital));
        //        assertEquals("DRUP Recip", HospitalModel.determinePIEmail("DRUP01010001", hospital));

        // Revert to CPCT PI with '*' for DRUP PI & recipients
      //  HospitalData hospital2 = hospitalModel.hospitalPerId("02");
//        assertNotNull(hospital2);
        // assertEquals("CPCT Recip2", HospitalModel.determinePIEmail("CPCT02010001", hospital2));
    }

    @Test
    public void extractHospitalName() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        //   assertEquals("N/A", hospitalModel.queryHospitalDataForSample("WIDE01010001").hospitalName());
        //    assertEquals("N/A", hospitalModel.queryHospitalDataForSample("CORE01010001").hospitalName());
    }

    @Test
    public void extractPIName() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        //   assertEquals("N/A", hospitalModel.queryHospitalDataForSample("WIDE01010001").analyseRequestName());
    }

    @Test
    public void extractPIEmail() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        //   assertEquals("N/A", hospitalModel.queryHospitalDataForSample("WIDE01010001").analyseRequestEmail());
    }

    @Test
    public void canReadHospitalNameAndAddress() {
        HospitalModel hospitalModel = buildTestHospitalModel();
     //   HospitalData hospital = hospitalModel.hospitalPerId("01");

      //  assertNotNull(hospital);
      //  assertEquals("ExtHosp1", hospital.externalHospitalName());
      //  assertEquals("Zip", hospital.addressZip());
      //  assertEquals("City", hospital.addressCity());
    }

    @Test
    public void canLookupAddresseeForSample() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        //  assertEquals("N/A", hospitalModel.queryHospitalDataForSample("CPCT02010001T").hospitalAdres());
    }

    @Test
    public void canLookupAddressForCORESample() {
        HospitalModel hospitalModel = buildTestHospitalModel();
        //  assertEquals("N/A", hospitalModel.queryHospitalDataForSample("CORE18001224T").hospitalAdres());
    }

    @NotNull
    private static HospitalModel buildTestHospitalModel() {
        Map<String, HospitalSampleMapping> hospitalPerIdManual = Maps.newHashMap();

        hospitalPerIdManual.put("CORE18001224T", ImmutableHospitalSampleMapping.of("HOSP1"));

        return ImmutableHospitalModel.of(
                hospitalPerIdManual,
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap());
    }
}
