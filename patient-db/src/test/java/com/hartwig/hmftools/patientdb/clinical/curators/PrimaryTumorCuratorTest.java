package com.hartwig.hmftools.patientdb.clinical.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;

import org.junit.Test;

public class PrimaryTumorCuratorTest {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canDetermineUnusedTerms() {
        PrimaryTumorCurator curator = CuratorTestFactory.primaryTumorCurator();
        assertEquals(5, curator.unusedSearchTerms().size());

        curator.search("patient", "desmoïd tumor");
        assertEquals(4, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // See DEV-275
        PrimaryTumorCurator curator = CuratorTestFactory.primaryTumorCurator();
        String desmoidTumor = "desmoïd tumor";
        CuratedPrimaryTumor primaryTumor = curator.search("patient", desmoidTumor);

        assertEquals("Bone/Soft tissue", primaryTumor.location());
    }

    @Test
    public void canOverrideTumorLocation() {
        PrimaryTumorCurator curator = CuratorTestFactory.primaryTumorCurator();
        CuratedPrimaryTumor primaryTumor = curator.search("PT1", "Does not curate");

        assertEquals("Adrenal gland", primaryTumor.location());
    }

    @Test
    public void canResolveDoidNodes() throws IOException {
        List<DoidNode> doidNodes = DiseaseOntology.readDoidOwlEntryFromDoidJson(DOID_FILE_JSON).nodes();
        assertEquals(Lists.newArrayList(doidNodes.get(0)), PrimaryTumorCurator.resolveDoidNodes(doidNodes, Lists.newArrayList("8718")));
    }

    @Test
    public void canCurateSearchTermWithChar34() {
        String searchTerm = "Non-small cell carcinoma NOS (mostly resembling lung carcinoma): working diagnosis \"lung carcinoma\"";
        PrimaryTumorCurator curator = CuratorTestFactory.primaryTumorCurator();
        CuratedPrimaryTumor primaryTumor = curator.search("patient", searchTerm);

        String location = primaryTumor.location();
        assertNotNull(location);
        assertEquals("lung", location.toLowerCase());
    }
}