package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;

import org.junit.Test;

public class TumorLocationCuratorV2Test {

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        assertEquals(3, curator.unusedSearchTerms().size());

        curator.search("Morbus Kahler");
        assertEquals(2, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // See DEV-275
        TumorLocationCuratorV2 curator = TestCuratorFactory.tumorLocationV2Curator();
        String desmoidTumor = "desmoïd tumor";
        CuratedTumorLocationV2 tumorLocation = curator.search(desmoidTumor);

        assertEquals("Bone/Soft tissue", tumorLocation.primaryTumorLocation());
    }
}