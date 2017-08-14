package com.hartwig.hmftools.common.centra;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class CentraTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("centra").getPath();
    private static final String TEST_FILE = BASE_RESOURCE_DIR + File.separator + "centra.csv";

    @Test
    public void canReadCPCTRecipients() throws IOException, EmptyFileException {
        final Map<String, String> recipientsPerCentra = Centra.readCPCTRecipientsFromCSV(TEST_FILE);
        assertEquals("my@email.com; my2@email.com", recipientsPerCentra.get("01"));
    }

    @Test
    public void canReadDRUPRecipients() throws IOException, EmptyFileException {
        final Map<String, String> recipientsPerCentra = Centra.readDRUPRecipientsFromCSV(TEST_FILE);
        assertEquals("my3@email.com", recipientsPerCentra.get("01"));
    }
}
