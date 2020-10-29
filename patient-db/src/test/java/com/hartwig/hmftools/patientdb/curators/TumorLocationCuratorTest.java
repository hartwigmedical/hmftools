package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;

import org.junit.Test;

public class TumorLocationCuratorTest {

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        assertEquals(5, curator.unusedSearchTerms().size());

        curator.search("Breast cancer");
        assertEquals(5, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // See DEV-275
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        String desmoidTumor = "desmoïd tumor";
        CuratedTumorLocationV2 tumorLocation = curator.search(desmoidTumor);

        String location = tumorLocation.primaryTumorLocation();
        assertNotNull(location);
        assertEquals("bone/soft tissue", location.toLowerCase());
    }

    @Test
    public void canCurateSearchTermWithChar34() {
        String searchTerm = "Non-small cell carcinoma NOS (mostly resembling lung carcinoma): working diagnosis \"lung carcinoma\"";
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        CuratedTumorLocationV2 tumorLocation = curator.search(searchTerm);

        String location = tumorLocation.primaryTumorLocation();
        assertNotNull(location);
        assertEquals("lung", location.toLowerCase());
    }
}