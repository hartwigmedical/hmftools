package com.hartwig.hmftools.common.centra;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class CentraTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("centra").getPath();
    private static final String TEST_FILE = BASE_RESOURCE_DIR + File.separator + "centra.csv";

    @Test
    public void canReadCPCTRecipients() throws IOException, EmptyFileException {
        final CentraModel centraModel = Centra.readFromCSV(TEST_FILE);
        assertEquals("my@email.com; my2@email.com", centraModel.getCpctRecipients("01"));
    }

    @Test
    public void canReadDRUPRecipients() throws IOException, EmptyFileException {
        final CentraModel centraModel = Centra.readFromCSV(TEST_FILE);
        assertEquals("my3@email.com", centraModel.getDrupRecipients("01"));
    }

    @Test
    public void canReadAddress() throws IOException, EmptyFileException {
        final CentraModel centraModel = Centra.readFromCSV(TEST_FILE);
        final CentraData centra = centraModel.centraPerId("01");
        assertNotNull(centra);
        assertEquals("Address-AVL", centra.addressName());
        assertEquals("1000 AB", centra.addressZip());
        assertEquals("AMSTERDAM", centra.addressCity());
    }
}
