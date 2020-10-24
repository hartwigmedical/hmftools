package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DiseaseOntologyTest {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canLoadDoidJsonFile() throws IOException {
        List<DoidEntry> doidEntries = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON);
        assertEquals(2, doidEntries.size());

        DoidEntry doidEntry1 = doidEntries.get(0);
        assertEquals(doidEntry1.url(), "http://purl.obolibrary.org/obo/DOID_8718");
        assertEquals(doidEntry1.doid(), "8718");
        assertEquals(doidEntry1.doidTerm(), "obsolete carcinoma in situ of respiratory system");

        DoidEntry doidEntry2 = doidEntries.get(1);
        assertEquals(doidEntry2.url(), "http://purl.obolibrary.org/obo/DOID_8717");
        assertEquals(doidEntry2.doid(), "8717");
        assertEquals(doidEntry2.doidTerm(), "decubitus ulcer");
    }
}