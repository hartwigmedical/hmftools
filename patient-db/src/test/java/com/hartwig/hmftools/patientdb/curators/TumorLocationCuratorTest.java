package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;

import org.junit.Test;

public class TumorLocationCuratorTest {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCurator curator = TestCuratorFactory.tumorLocationCurator();
        assertEquals(5, curator.unusedSearchTerms().size());

        curator.search("Morbus Kahler");
        assertEquals(4, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // See DEV-275
        TumorLocationCurator curator = TestCuratorFactory.tumorLocationCurator();
        String desmoidTumor = "desmoïd tumor";
        CuratedTumorLocation tumorLocation = curator.search(desmoidTumor);

        assertEquals("Bone/Soft tissue", tumorLocation.primaryTumorLocation());
    }

    @Test
    public void canResolveDoidNodes() throws IOException {
        List<DoidNode> doidNodes = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON).nodes();
        List<String> doids = Lists.newArrayList();
        doids.add("8718");
        assertEquals(Lists.newArrayList(doidNodes.get(0)), TumorLocationCurator.resolveDoidNodes(doidNodes, doids));
    }

    @Test
    public void canCurateSearchTermWithChar34() {
        String searchTerm = "Non-small cell carcinoma NOS (mostly resembling lung carcinoma): working diagnosis \"lung carcinoma\"";
        TumorLocationCurator curator = TestCuratorFactory.tumorLocationCurator();
        CuratedTumorLocation tumorLocation = curator.search(searchTerm);

        String location = tumorLocation.primaryTumorLocation();
        assertNotNull(location);
        assertEquals("lung", location.toLowerCase());
    }
}