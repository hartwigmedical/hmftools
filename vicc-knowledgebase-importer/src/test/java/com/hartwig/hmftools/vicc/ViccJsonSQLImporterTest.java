package com.hartwig.hmftools.vicc;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.junit.Ignore;
import org.junit.Test;

public class ViccJsonSQLImporterTest {

    private static final String VICC_JSON_FILE = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

    @Test
    @Ignore
    public void readViccJson() throws IOException {
        // This function exists just for fast local testing.
        List<ViccEntry> entries = ViccJsonReader.readAll(VICC_JSON_FILE);
        assertNotNull(entries);
    }
}
