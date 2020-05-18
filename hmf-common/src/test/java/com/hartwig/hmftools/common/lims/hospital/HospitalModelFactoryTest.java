package com.hartwig.hmftools.common.lims.hospital;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canCreateFromLimsDirectoryWithoutWarnings() throws IOException {
        HospitalModel hospitalModel = HospitalModelFactory.fromLimsDirectory(LIMS_DIRECTORY);
        assertNotNull(hospitalModel);

        assertTrue(HospitalModelFactory.validateModelIntegrity(hospitalModel));
    }

    @Test
    public void canValidateIntegrityCorrectly() {
        HospitalModel emptyModel = ImmutableHospitalModel.builder().build();
        assertTrue(HospitalModelFactory.validateModelIntegrity(emptyModel));

        HospitalModel noAddressModel = ImmutableHospitalModel.builder()
                .putHospitalPersonsCPCT("HOSP",
                        ImmutableHospitalPersons.builder()
                                .hospitalPI("Test")
                                .requesterName("Test")
                                .requesterEmail("Test")
                                .build())
                .build();

        assertFalse(HospitalModelFactory.validateModelIntegrity(noAddressModel));
    }

    @Test
    public void canReadHospitalAddressTsv() throws IOException {
        Map<String, HospitalAddress> hospitalAddress =
                HospitalModelFactory.readFromHospitalAddress(LIMS_DIRECTORY + File.separator + "hospital_address.tsv");

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
                HospitalModelFactory.readFromHospitalPersons(LIMS_DIRECTORY + File.separator + "hospital_cpct.tsv", 2);

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
                HospitalModelFactory.readFromHospitalPersons(LIMS_DIRECTORY + File.separator + "hospital_drup.tsv", 2);
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
    public void canReadHospitalContactWIDE() throws IOException {
        Map<String, HospitalPersons> hospitalContactWIDE =
                HospitalModelFactory.readFromHospitalPersons(LIMS_DIRECTORY + File.separator + "hospital_wide.tsv", 4);
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
    public void canReadSampleToHospitalMapping() throws IOException {
        Map<String, String> sampleHospitalMapping =
                HospitalModelFactory.readFromSampleToHospitalMapping(LIMS_DIRECTORY + File.separator + "sample_hospital_mapping.tsv");

        assertEquals(2, sampleHospitalMapping.size());

        assertEquals("01", sampleHospitalMapping.get("CORE18001224T"));
        assertEquals("02", sampleHospitalMapping.get("CORE18002000T"));
    }
}
