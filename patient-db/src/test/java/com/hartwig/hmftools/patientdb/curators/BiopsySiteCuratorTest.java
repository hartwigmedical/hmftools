package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.FileInputStream;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientdb.data.CuratedBiopsyType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsySiteCuratorTest {

    private static final String BIOPSY_SITE_MAPPING_CSV = Resources.getResource("test_biopsy_site_mapping.csv").getPath();

    @Test
    public void canCreateFromProductionResource() throws IOException {
        assertNotNull(BiopsySiteCurator.fromProductionResource());
    }

    @Test
    public void canCurateOnTestResource() {
        BiopsySiteCurator curator = createBiopsySiteCurator();

        CuratedBiopsyType knownCuratedType = curator.search("Breast", "Primary");
        assertEquals("Breast", knownCuratedType.type());

        CuratedBiopsyType cannotCurate = curator.search("This is", "Unknown");
        assertNull(cannotCurate.type());
    }

    @NotNull
    private static BiopsySiteCurator createBiopsySiteCurator() {
        try {
            return new BiopsySiteCurator(new FileInputStream(BIOPSY_SITE_MAPPING_CSV));
        } catch (IOException e) {
            throw new IllegalStateException("Could not create biopsy site curator!");
        }
    }
}