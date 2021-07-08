package com.hartwig.hmftools.common.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canReadHospitalAddressTsv() throws IOException {
        Map<String, HospitalAddress> hospitalAddress =
                HospitalModelFactory.readFromHospitalAddressTsv(LIMS_DIRECTORY + File.separator + "hospital_address.tsv");

        assertEquals(2, hospitalAddress.size());

        HospitalAddress address1 = hospitalAddress.get("01");
        assertEquals("Ext-HMF", address1.hospitalName());
        assertEquals("1000 AB", address1.hospitalZip());
        assertEquals("AMSTERDAM", address1.hospitalCity());

        HospitalAddress address2 = hospitalAddress.get("02");
        assertEquals("Ext-HMF", address2.hospitalName());
        assertEquals("1000 AB", address2.hospitalZip());
        assertEquals("AMSTERDAM", address2.hospitalCity());
    }

    @Test
    public void canReadHospitalContactCPCT() throws IOException {
        Map<String, HospitalPersons> hospitalContactCPCT =
                HospitalModelFactory.readFromHospitalPersonsTsv(LIMS_DIRECTORY + File.separator + "hospital_cpct.tsv", 2, "CPCT");

        assertEquals(2, hospitalContactCPCT.size());

        HospitalPersons cpct1 = hospitalContactCPCT.get("01");
        assertEquals("Someone", cpct1.hospitalPI());
        assertNull(cpct1.requesterName());
        assertNull(cpct1.requesterEmail());

        HospitalPersons cpct2 = hospitalContactCPCT.get("02");
        assertEquals("Someone", cpct2.hospitalPI());
        assertNull(cpct2.requesterName());
        assertNull(cpct2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactDRUP() throws IOException {
        Map<String, HospitalPersons> hospitalContactDRUP =
                HospitalModelFactory.readFromHospitalPersonsTsv(LIMS_DIRECTORY + File.separator + "hospital_drup.tsv", 2, "DRUP");
        assertEquals(2, hospitalContactDRUP.size());

        HospitalPersons drup1 = hospitalContactDRUP.get("01");
        assertEquals("Someone", drup1.hospitalPI());
        assertNull(drup1.requesterName());
        assertNull(drup1.requesterEmail());

        HospitalPersons drup2 = hospitalContactDRUP.get("02");
        assertEquals("Someone", drup2.hospitalPI());
        assertNull(drup2.requesterName());
        assertNull(drup2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactACTIN() throws IOException {
        Map<String, HospitalPersons> hospitalContactACTIN =
                HospitalModelFactory.readFromHospitalPersonsTsv(LIMS_DIRECTORY + File.separator + "hospital_actin.tsv", 2, "ACTIN");
        assertEquals(1, hospitalContactACTIN.size());

        HospitalPersons actin2 = hospitalContactACTIN.get("02");
        assertEquals("Someone", actin2.hospitalPI());
        assertNull(actin2.requesterName());
        assertNull(actin2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactWIDE() throws IOException {
        Map<String, HospitalPersons> hospitalContactWIDE =
                HospitalModelFactory.readFromHospitalPersonsTsv(LIMS_DIRECTORY + File.separator + "hospital_wide.tsv", 4, "WIDE");
        assertEquals(2, hospitalContactWIDE.size());

        HospitalPersons wide1 = hospitalContactWIDE.get("01");
        assertEquals("Someone", wide1.hospitalPI());
        assertEquals("Someone1", wide1.requesterName());
        assertEquals("my@email.com", wide1.requesterEmail());

        HospitalPersons wide2 = hospitalContactWIDE.get("02");
        assertEquals("Someone", wide2.hospitalPI());
        assertEquals("Someone1", wide2.requesterName());
        assertEquals("my@email.com", wide2.requesterEmail());
    }

    @Test
    public void canReadHospitalContactCOREDB() throws IOException {
        Map<String, HospitalPersons> hospitalContactCOREDB =
                HospitalModelFactory.readFromHospitalPersonsTsv(LIMS_DIRECTORY + File.separator + "hospital_coredb.tsv", 4, "COREDB");
        assertEquals(1, hospitalContactCOREDB.size());

        HospitalPersons coredb01 = hospitalContactCOREDB.get("01");
        assertEquals("Someone1", coredb01.hospitalPI());
        assertEquals("Someone1", coredb01.requesterName());
        assertEquals("my@email.com", coredb01.requesterEmail());
    }

    @Test
    public void canReadSampleToHospitalMapping() throws IOException {
        Map<String, String> sampleHospitalMapping =
                HospitalModelFactory.readFromSampleToHospitalMappingTsv(LIMS_DIRECTORY + File.separator + "sample_hospital_mapping.tsv");

        assertEquals(2, sampleHospitalMapping.size());

        assertEquals("01", sampleHospitalMapping.get("CORE18001224T"));
        assertEquals("02", sampleHospitalMapping.get("CORE18002000T"));
    }
}
