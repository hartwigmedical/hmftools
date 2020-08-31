package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodFile;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverGeneLikelihood;

import org.jetbrains.annotations.NotNull;

public class DriverGenePanelFactory {

    @Deprecated
    public DriverGenePanel create() {
        return create(builtIn());
    }

    @NotNull
    private static List<DriverGene> builtIn() {
        final InputStream inputStream = DriverGenePanel.class.getResourceAsStream("/drivercatalog/DriverGenePanel.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("gene"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @NotNull
    public static DriverGenePanel empty() {
        return create(Collections.emptyList());
    }

    @NotNull
    public static DriverGenePanel buildFromTsv(@NotNull String driverGenePanelFile) throws IOException {
        return create(DriverGeneFile.read(driverGenePanelFile));
    }

    public static DriverGenePanel create(final List<DriverGene> genes) {
        return create(DriverGenePanelAssembly.HG19, genes);
    }

    @NotNull
    static DriverGenePanel create(@NotNull final DriverGenePanelAssembly assembly, final List<DriverGene> genes) {
        final DndsGeneName dndsGeneName = new DndsGeneName(assembly);
        final Map<String, String> dndsTsGenes = Maps.newHashMap();
        final Map<String, String> dndsOncoGenes = Maps.newHashMap();

        for (DriverGene gene : genes) {
            if (gene.reportVariant()) {
                boolean isValidGene = dndsGeneName.isValid(gene);
                if (!isValidGene) {
                    throw new IllegalArgumentException(
                            "dNdS not available for " + gene.likelihoodType() + " gene " + gene.gene() + " in driver gene panel");
                }

                final String dndsName = dndsGeneName.dndsGeneName(gene);
                if (gene.likelihoodType().equals(DriverCategory.TSG)) {
                    dndsTsGenes.put(dndsName, gene.gene());
                } else {
                    dndsOncoGenes.put(dndsName, gene.gene());
                }
            }
        }

        Map<String, DndsDriverGeneLikelihood> tsgLikelihood = tsgLikelihood().stream()
                .filter(x -> dndsTsGenes.containsKey(x.gene()))
                .map(x -> ImmutableDndsDriverGeneLikelihood.builder().from(x).gene(dndsTsGenes.get(x.gene())).build())
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));
        Map<String, DndsDriverGeneLikelihood> oncoLikelihood = oncoLikelihood().stream()
                .filter(x -> dndsOncoGenes.containsKey(x.gene()))
                .map(x -> ImmutableDndsDriverGeneLikelihood.builder().from(x).gene(dndsOncoGenes.get(x.gene())).build())
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        Map<String, String> deletionBandMap = genes.stream()
                .filter(x -> x.reportDeletion() && !x.deletionBand().isEmpty())
                .collect(Collectors.toMap(DriverGene::gene, DriverGene::deletionBand));

        return ImmutableDriverGenePanel.builder()
                .driverGenes(genes)
                .tsgLikelihood(tsgLikelihood)
                .oncoLikelihood(oncoLikelihood)
                .deletionBandMap(deletionBandMap)
                .build();
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> tsgLikelihood() {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/DndsDriverLikelihoodTsg.tsv"));
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> oncoLikelihood() {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/DndsDriverLikelihoodOnco.tsv"));
    }

    @NotNull
    private static InputStream resource(@NotNull final String resource) {
        return DriverGenePanelFactory.class.getResourceAsStream(resource);
    }
}
