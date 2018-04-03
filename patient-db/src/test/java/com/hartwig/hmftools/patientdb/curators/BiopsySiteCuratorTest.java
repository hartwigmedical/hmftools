package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import com.hartwig.hmftools.patientdb.data.CuratedBiopsyType;

import org.junit.Test;

public class BiopsySiteCuratorTest {

    @Test
    public void canCreateFromProductionResource() throws IOException {
        assertNotNull(BiopsySiteCurator.fromProductionResource());
    }

    @Test
    public void canCurateOnTestResource() {
        BiopsySiteCurator curator = TestCuratorFactory.biopsySiteCurator();

        CuratedBiopsyType knownCuratedType = curator.search("Breast", "HER2 Positive", "Primary", "X");
        assertEquals("Breast", knownCuratedType.type());

        CuratedBiopsyType fallbackCuratedType = curator.search("Breast", "HER2 Positive", "Primary", null);
        assertEquals("Breast", fallbackCuratedType.type());

        CuratedBiopsyType cannotCurate = curator.search("This", "Is", "Unknown", "X");
        assertNull(cannotCurate.type());

        CuratedBiopsyType cannotCurate2 = curator.search(null, null, null, null);
        assertNull(cannotCurate2.type());
    }
}