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
        List<Doid> doids = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON);
        assertEquals(2, doids.size());

        Doid doid1 = doids.get(0);
        assertEquals(doid1.id(), "http://purl.obolibrary.org/obo/DOID_8718");
        assertEquals(doid1.doid(), "8718");
        assertEquals(doid1.doidTerm(), "obsolete carcinoma in situ of respiratory system");

        Doid doid2 = doids.get(1);
        assertEquals(doid2.id(), "http://purl.obolibrary.org/obo/DOID_8717");
        assertEquals(doid2.doid(), "8717");
        assertEquals(doid2.doidTerm(), "decubitus ulcer");
    }
}