package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.FileInputStream;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientdb.data.CuratedCancerType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorLocationCuratorTest {

    private static final String TUMOR_LOCATION_MAPPING_CSV = Resources.getResource("test_tumor_location_mapping.csv").getPath();

    @Test
    public void canCreateFromProductionResource() throws IOException {
        assertNotNull(TumorLocationCurator.fromProductionResource());
    }

    @Test
    public void canDetermineUnusedTerms() {
        TumorLocationCurator curator = createTumorLocationCurator();
        assertEquals(7, curator.unusedSearchTerms().size());

        curator.search("Breast cancer");
        assertEquals(6, curator.unusedSearchTerms().size());
    }

    @Test
    public void canCurateDesmoidTumor() {
        // KODU: See DEV-275
        TumorLocationCurator curator = createTumorLocationCurator();
        String desmoidTumor = "desmoïd tumor";
        CuratedCancerType cancerType = curator.search(desmoidTumor);

        String category = cancerType.category();
        assertNotNull(category);
        assertEquals("sarcoma", category.toLowerCase());
    }

    @NotNull
    private static TumorLocationCurator createTumorLocationCurator() {
        try {
            return new TumorLocationCurator(new FileInputStream(TUMOR_LOCATION_MAPPING_CSV));
        } catch (IOException e) {
            throw new IllegalStateException("Could not create tumor location curator!");
        }
    }
}