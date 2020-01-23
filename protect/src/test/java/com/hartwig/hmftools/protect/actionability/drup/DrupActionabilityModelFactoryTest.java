package com.hartwig.hmftools.protect.actionability.drup;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class DrupActionabilityModelFactoryTest {

    private static final String DRUP_GENES_CSV = Resources.getResource("actionability/drup_genes.csv").getPath();

    @Test
    public void canReadTestDrupCsv() throws IOException {
        DrupActionabilityModel drupActionabilityModel = DrupActionabilityModelFactory.buildFromCsv(DRUP_GENES_CSV);

        assertEquals(2, drupActionabilityModel.actionableGenes().size());
    }
}