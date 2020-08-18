package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;

import org.jetbrains.annotations.NotNull;

public class DriverGenePanel {

    private final List<DriverGene> genes;
    private final Map<String, DndsDriverGeneLikelihood> tsgLikelihood;
    private final Map<String, DndsDriverImpactLikelihood> oncoLikelihood;
    private final Set<String> amplificationTargets;
    private final Set<String> deletionTargets;
    private final Map<String, String> deletionBandMap;


    public DriverGenePanel(final List<DriverGene> genes) {
        this.genes = genes;
        final Set<String> tsGenes = genes.stream()
                .filter(x -> x.likelihoodType().equals(DriverLikelihoodType.TSG))
                .map(DriverGene::gene)
                .collect(Collectors.toSet());

        final Set<String> oncoGenes = genes.stream()
                .filter(x -> x.likelihoodType().equals(DriverLikelihoodType.ONCO))
                .map(DriverGene::gene)
                .collect(Collectors.toSet());

        tsgLikelihood = DndsDriverGeneLikelihoodSupplier.tsgLikelihood(tsGenes);
        oncoLikelihood = DndsDriverGeneLikelihoodSupplier.oncoLikelihood(oncoGenes);
        amplificationTargets = genes.stream().filter(DriverGene::reportAmplification).map(DriverGene::gene).collect(Collectors.toSet());
        deletionTargets = genes.stream().filter(DriverGene::reportDeletionAndDisruption).map(DriverGene::gene).collect(Collectors.toSet());
        deletionBandMap = genes.stream()
                .filter(x -> x.reportDeletionAndDisruption() && !x.deletionBand().isEmpty())
                .collect(Collectors.toMap(DriverGene::gene, DriverGene::deletionBand));
    }

    @NotNull
    public Set<String> oncoGenes() {
        return oncoLikelihood.keySet();
    }

    @NotNull
    public Set<String> tsGenes() {
        return tsgLikelihood.keySet();
    }

    @NotNull
    public Map<String, DndsDriverGeneLikelihood> tsgLikelihood() {
        return tsgLikelihood;
    }

    @NotNull
    public Map<String, DndsDriverImpactLikelihood> oncoLikelihood() {
        return oncoLikelihood;
    }

    public Set<String> amplificationTargets() {
        return amplificationTargets;
    }

    public Set<String> deletionTargets() {
        return deletionTargets;
    }

    public Map<String, String> deletionBandMap() {
        return deletionBandMap;
    }
}
