package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;

import org.junit.Test;

public class TumorLocationCuratorV2Test {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        assertEquals(5, curator.unusedSearchTerms().size());

        curator.search("Morbus Kahler");
        assertEquals(4, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // See DEV-275
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        String desmoidTumor = "desmoïd tumor";
        CuratedTumorLocationV2 tumorLocation = curator.search(desmoidTumor);

        assertEquals("Bone/Soft tissue", tumorLocation.primaryTumorLocation());
    }

    @Test
    public void canResolveDoidNodes() throws IOException {
        List<DoidNode> doidNodes = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON).nodes();
        List<String> doids = Lists.newArrayList();
        doids.add("8718");
        assertEquals(Lists.newArrayList(doidNodes.get(0)), TumorLocationCuratorV2.resolveDoidNodes(doidNodes, doids));
    }
}