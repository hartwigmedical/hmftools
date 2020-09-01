package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.InputStream;
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

    @NotNull
    public static DriverGenePanel empty() {
        return create(DriverGenePanelAssembly.HG19, Collections.emptyList());
    }

    @NotNull
    public static DriverGenePanel create(@NotNull final DriverGenePanelAssembly assembly, final List<DriverGene> genes) {
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

        return ImmutableDriverGenePanel.builder()
                .driverGenes(genes)
                .tsgLikelihood(tsgLikelihood)
                .oncoLikelihood(oncoLikelihood)
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
