package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canReadFromHospitalDirectory() throws IOException {
        HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(LIMS_DIRECTORY);
        assertNotNull(hospitalModel);
    }

    @Test
    public void canReadHospitalAdress() throws IOException {
        Map<String, HospitalAdress> hospitalAdress =
                HospitalModelFactory.readFromHospitalAdress(LIMS_DIRECTORY + File.separator + "hospital_adress.csv");

        assertEquals(2, hospitalAdress.size());

        HospitalAdress adress1 = hospitalAdress.get("01");
        assertEquals("01", adress1.hospitalId());
        assertEquals("Ext-HMF", adress1.hospitalName());
        assertEquals("1000 AB", adress1.hospitalZip());
        assertEquals("AMSTERDAM", adress1.hospitalCity());

        HospitalAdress adress2 = hospitalAdress.get("02");
        assertEquals("02", adress2.hospitalId());
        assertEquals("Ext-HMF", adress2.hospitalName());
        assertEquals("1000 AB", adress2.hospitalZip());
        assertEquals("AMSTERDAM", adress2.hospitalCity());
    }

    @Test
    @Ignore
    // TODO Will be fixed once requester columns have been removed.
    public void canReadHospitalCPCT() throws IOException {
        Map<String, HospitalData> hospitalDataCPCT =
                HospitalModelFactory.readFromHospitalDataCPCT(LIMS_DIRECTORY + File.separator + "hospital_cpct.csv");

        assertEquals(2, hospitalDataCPCT.size());

        HospitalData cpct1 = hospitalDataCPCT.get("01");
        assertEquals("01", cpct1.hospitalId());
        assertEquals("Someone", cpct1.hospitalPI());
        assertEquals("", cpct1.requestName());
        assertEquals("", cpct1.requestEmail());

        HospitalData cpct2 = hospitalDataCPCT.get("01");
        assertEquals("01", cpct2.hospitalId());
        assertEquals("Someone", cpct2.hospitalPI());
        assertEquals("", cpct2.requestName());
        assertEquals("", cpct2.requestEmail());

    }

    @Test
    @Ignore
    // TODO Will be fixed once requester columns have been removed.
    public void canReadHospitalDRUP() throws IOException {
        Map<String, HospitalData> hospitalDataDRUP =
                HospitalModelFactory.readFromHospitalDataDRUP(LIMS_DIRECTORY + File.separator + "hospital_drup.csv");
        assertEquals(2, hospitalDataDRUP.size());

        HospitalData drup1 = hospitalDataDRUP.get("01");
        assertEquals("01", drup1.hospitalId());
        assertEquals("Someone", drup1.hospitalPI());
        assertEquals("", drup1.requestName());
        assertEquals("", drup1.requestEmail());

        HospitalData drup2 = hospitalDataDRUP.get("01");
        assertEquals("01", drup2.hospitalId());
        assertEquals("Someone", drup2.hospitalPI());
        assertEquals("", drup2.requestName());
        assertEquals("", drup2.requestEmail());

    }

    @Test
    public void canReadHospitalWIDE() throws IOException {
        Map<String, HospitalData> hospitalDataWIDE =
                HospitalModelFactory.readFromHospitalDataWIDE(LIMS_DIRECTORY + File.separator + "hospital_wide.csv");
        assertEquals(2, hospitalDataWIDE.size());

        HospitalData wide1 = hospitalDataWIDE.get("01");
        assertEquals("01", wide1.hospitalId());
        assertEquals("Someone", wide1.hospitalPI());
        assertEquals("Someone1", wide1.requestName());
        assertEquals("my@email.com", wide1.requestEmail());

        HospitalData wide2 = hospitalDataWIDE.get("02");
        assertEquals("02", wide2.hospitalId());
        assertEquals("Someone", wide2.hospitalPI());
        assertEquals("Someone1", wide2.requestName());
        assertEquals("my@email.com", wide2.requestEmail());
    }

    @Test
    public void canReadSampleHospitalMapping() throws IOException {
        Map<String, HospitalSampleMapping> sampleHospitalMapping =
                HospitalModelFactory.readFromSampleHospitalMapping(LIMS_DIRECTORY + File.separator + "sample_hospital_mapping.csv");

        assertEquals(1, sampleHospitalMapping.size());

        HospitalSampleMapping sampleMapping = sampleHospitalMapping.get("CORE18001224T");
        assertEquals("HOSP1", sampleMapping.internalHospitalName());
    }
}
