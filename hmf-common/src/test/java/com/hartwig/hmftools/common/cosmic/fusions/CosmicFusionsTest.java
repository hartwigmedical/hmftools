package com.hartwig.hmftools.common.cosmic.fusions;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CosmicFusionsTest {

    private static final String FUSION_EXAMPLE_FILE = Resources.getResource("cosmic").getPath() + File.separator + "FusionsExample.csv";

    @Test
    public void canReadFromCSV() throws IOException {
        CosmicFusionModel fusionModel = CosmicFusions.readFromCSV(FUSION_EXAMPLE_FILE);
        assertEquals(5, fusionModel.fusions().size());

        assertEquals(1, fusionModel.promiscuousFivePrime().size());
        assertEquals(0, fusionModel.promiscuousThreePrime().size());
    }
}
