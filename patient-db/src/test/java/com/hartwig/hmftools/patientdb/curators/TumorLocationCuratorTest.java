package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;

import org.junit.Test;

public class TumorLocationCuratorTest {

    @Test
    public void canCreateFromProductionResource() throws IOException {
        assertNotNull(TumorLocationCurator.fromProductionResource());
    }

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCurator curator = TestCuratorFactory.tumorLocationCurator();
        assertEquals(8, curator.unusedSearchTerms().size());

        curator.search("Breast cancer");
        assertEquals(7, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // KODU: See DEV-275
        TumorLocationCurator curator = TestCuratorFactory.tumorLocationCurator();
        String desmoidTumor = "desmoïd tumor";
        CuratedTumorLocation tumorLocation = curator.search(desmoidTumor);

        String location = tumorLocation.primaryTumorLocation();
        assertNotNull(location);
        assertEquals("sarcoma", location.toLowerCase());
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