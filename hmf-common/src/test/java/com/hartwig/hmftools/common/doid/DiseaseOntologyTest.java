package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DiseaseOntologyTest {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canExtractDoidFromUrl() {
        String url = "http://purl.obolibrary.org/obo/DOID_345";
        assertEquals("345", DiseaseOntology.extractDoid(url));
    }
    @Test
    public void canLoadDoidJsonFile() throws IOException {
        List<DoidNode> doidNodes = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON).doidNodes();
        assertEquals(2, doidNodes.size());

        DoidNode doidEntry1 = doidNodes.get(0);
        assertEquals(doidEntry1.url(), "http://purl.obolibrary.org/obo/DOID_8718");
        assertEquals(doidEntry1.doid(), "8718");
        assertEquals(doidEntry1.doidTerm(), "obsolete carcinoma in situ of respiratory system");

        DoidNode doidEntry2 = doidNodes.get(1);
        assertEquals(doidEntry2.url(), "http://purl.obolibrary.org/obo/DOID_8717");
        assertEquals(doidEntry2.doid(), "8717");
        assertEquals(doidEntry2.doidTerm(), "decubitus ulcer");
    }
}