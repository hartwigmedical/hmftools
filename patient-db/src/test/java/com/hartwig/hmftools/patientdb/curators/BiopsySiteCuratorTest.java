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

        CuratedBiopsyType knownCuratedType = curator.search("Breast", "HER2 Positive", "Primary");
        assertEquals("Breast", knownCuratedType.type());

        CuratedBiopsyType cannotCurate = curator.search("This", "Is", "Unknown");
        assertNull(cannotCurate.type());
    }
}