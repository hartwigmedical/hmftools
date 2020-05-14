package com.hartwig.hmftools.common.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canCreateFromLimsDirectory() throws IOException {
        HospitalModel hospitalModel = HospitalModelFactory.fromLimsDirectory(LIMS_DIRECTORY);
        assertNotNull(hospitalModel);
    }

    @Test
    public void canReadHospitalAddressTsv() throws IOException {
        Map<String, HospitalAddress> hospitalAddress =
                HospitalModelFactory.readFromHospitalAddress(LIMS_DIRECTORY + File.separator + "hospital_address.tsv");

        assertEquals(2, hospitalAddress.size());

        HospitalAddress address1 = hospitalAddress.get("01");
        assertEquals("01", address1.hospitalId());
        assertEquals("Ext-HMF", address1.hospitalName());
        assertEquals("1000 AB", address1.hospitalZip());
        assertEquals("AMSTERDAM", address1.hospitalCity());

        HospitalAddress address2 = hospitalAddress.get("02");
        assertEquals("02", address2.hospitalId());
        assertEquals("Ext-HMF", address2.hospitalName());
        assertEquals("1000 AB", address2.hospitalZip());
        assertEquals("AMSTERDAM", address2.hospitalCity());
    }

    @Test
    public void canReadHospitalContactCPCT() throws IOException {
        Map<String, HospitalContact> hospitalContactCPCT =
                HospitalModelFactory.readFromHospitalContact(LIMS_DIRECTORY + File.separator + "hospital_cpct.tsv", 2);

        assertEquals(2, hospitalContactCPCT.size());

        HospitalContact cpct1 = hospitalContactCPCT.get("01");
        assertEquals("01", cpct1.hospitalId());
        assertEquals("Someone", cpct1.hospitalPI());
        assertNull(cpct1.requesterName());
        assertNull(cpct1.requesterEmail());

        HospitalContact cpct2 = hospitalContactCPCT.get("02");
        assertEquals("02", cpct2.hospitalId());
        assertEquals("Someone", cpct2.hospitalPI());
        assertNull(cpct2.requesterName());
        assertNull(cpct2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactDRUP() throws IOException {
        Map<String, HospitalContact> hospitalContactDRUP =
                HospitalModelFactory.readFromHospitalContact(LIMS_DIRECTORY + File.separator + "hospital_drup.tsv", 2);
        assertEquals(2, hospitalContactDRUP.size());

        HospitalContact drup1 = hospitalContactDRUP.get("01");
        assertEquals("01", drup1.hospitalId());
        assertEquals("Someone", drup1.hospitalPI());
        assertNull(drup1.requesterName());
        assertNull(drup1.requesterEmail());

        HospitalContact drup2 = hospitalContactDRUP.get("02");
        assertEquals("02", drup2.hospitalId());
        assertEquals("Someone", drup2.hospitalPI());
        assertNull(drup2.requesterName());
        assertNull(drup2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactWIDE() throws IOException {
        Map<String, HospitalContact> hospitalContactWIDE =
                HospitalModelFactory.readFromHospitalContact(LIMS_DIRECTORY + File.separator + "hospital_wide.tsv", 4);
        assertEquals(2, hospitalContactWIDE.size());

        HospitalContact wide1 = hospitalContactWIDE.get("01");
        assertEquals("01", wide1.hospitalId());
        assertEquals("Someone", wide1.hospitalPI());
        assertEquals("Someone1", wide1.requesterName());
        assertEquals("my@email.com", wide1.requesterEmail());

        HospitalContact wide2 = hospitalContactWIDE.get("02");
        assertEquals("02", wide2.hospitalId());
        assertEquals("Someone", wide2.hospitalPI());
        assertEquals("Someone1", wide2.requesterName());
        assertEquals("my@email.com", wide2.requesterEmail());
    }

    @Test
    public void canReadSampleHospitalMapping() throws IOException {
        Map<String, HospitalSampleMapping> sampleHospitalMapping =
                HospitalModelFactory.readFromSampleHospitalMapping(LIMS_DIRECTORY + File.separator + "sample_hospital_mapping.tsv");

        assertEquals(1, sampleHospitalMapping.size());

        HospitalSampleMapping sampleMapping = sampleHospitalMapping.get("CORE18001224T");
        assertEquals("01", sampleMapping.hospitalId());
    }
}
