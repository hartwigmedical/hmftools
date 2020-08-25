package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class DriverGenePanelFactoryTest {

    @NotNull
    private static List<DriverGene> builtIn() {
        final InputStream inputStream = DriverGenePanelFactoryTest.class.getResourceAsStream("/drivercatalog/sample.driver.gene.panel.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("gene"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @NotNull
    public static DriverGenePanel testGenePanel() {
        return DriverGenePanelFactory.create(builtIn());
    }


}
