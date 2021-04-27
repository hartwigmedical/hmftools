package com.hartwig.hmftools.common.drivercatalog.panel;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DriverGenePanelFactoryTest {

    @NotNull
    public static DriverGenePanel testGenePanel() {
        return DriverGenePanelFactory.create(RefGenomeVersion.V37, builtIn());
    }

    @NotNull
    private static List<DriverGene> builtIn() {
        final InputStream inputStream = DriverGenePanelFactoryTest.class.getResourceAsStream("/drivercatalog/driver.gene.panel.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("gene"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @Test
    public void canReadOncoGenes() {
        DriverGenePanel genePanel = testGenePanel();
        DndsDriverImpactLikelihood missense = genePanel.oncoLikelihood().get("ABL1").missense();

        assertEquals(0.0014228, missense.driversPerSample(), 1e-7);
        assertEquals(4.436149558484221E-7, missense.passengersPerMutation(), 1e-8);
    }

    @Test
    public void canReadTSGGenes() {
        DriverGenePanel genePanel = testGenePanel();
        DndsDriverGeneLikelihood gene = genePanel.tsgLikelihood().get("ACVR1B");
        DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ACVR1B", gene.gene());
        assertEquals(0.003, missense.driversPerSample(), 0.001);
        assertEquals(2e-07, missense.passengersPerMutation(), 1e-7);
    }
}
