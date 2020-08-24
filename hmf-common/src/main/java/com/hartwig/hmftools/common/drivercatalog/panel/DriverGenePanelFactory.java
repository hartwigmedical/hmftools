package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodSupplier;

import org.jetbrains.annotations.NotNull;

public class DriverGenePanelFactory {

    @NotNull
    private static List<DriverGene> builtIn() {
        final InputStream inputStream = DriverGenePanel.class.getResourceAsStream("/drivercatalog/DriverGenePanel.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("gene"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @Deprecated
    public DriverGenePanel create() {
        return create(builtIn());
    }

    public DriverGenePanel empty() {
        return create(Collections.emptyList());
    }

    @NotNull
    public DriverGenePanel create(final List<DriverGene> genes) {

        final Set<String> tsGenes = genes.stream()
                .filter(x -> x.likelihoodType().equals(DriverCategory.TSG) && x.reportVariant())
                .map(DriverGene::gene)
                .collect(Collectors.toSet());

        final Set<String> oncoGenes = genes.stream()
                .filter(x -> x.likelihoodType().equals(DriverCategory.ONCO) && x.reportVariant())
                .map(DriverGene::gene)
                .collect(Collectors.toSet());

        Map<String, DndsDriverGeneLikelihood> tsgLikelihood = DndsDriverGeneLikelihoodSupplier.tsgLikelihood(tsGenes);
        Map<String, DndsDriverGeneLikelihood> oncoLikelihood = DndsDriverGeneLikelihoodSupplier.oncoLikelihood(oncoGenes);
        Set<String> amplificationTargets =
                genes.stream().filter(DriverGene::reportAmplification).map(DriverGene::gene).collect(Collectors.toSet());
        Set<String> deletionTargets =
                genes.stream().filter(DriverGene::reportDeletionAndDisruption).map(DriverGene::gene).collect(Collectors.toSet());
        Map<String, String> deletionBandMap = genes.stream()
                .filter(x -> x.reportDeletionAndDisruption() && !x.deletionBand().isEmpty())
                .collect(Collectors.toMap(DriverGene::gene, DriverGene::deletionBand));

        return ImmutableDriverGenePanel.builder()
                .driverGenes(genes)
                .tsgLikelihood(tsgLikelihood)
                .oncoLikelihood(oncoLikelihood)
                .amplificationTargets(amplificationTargets)
                .deletionTargets(deletionTargets)
                .deletionBandMap(deletionBandMap)
                .build();
    }

}
