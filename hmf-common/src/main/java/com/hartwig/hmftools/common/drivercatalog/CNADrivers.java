package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

public class CNADrivers {

    private static final double MIN_COPY_NUMBER_RELATIVE_INCREASE = 3;
    private static final double MAX_COPY_NUMBER_DEL = 0.5;

    @NotNull
    private final Set<String> oncoGenes;
    @NotNull
    private final Set<String> tsGenes;
    @NotNull
    private final Set<String> amplificationTargets;
    @NotNull
    private final Set<String> deletionTargets;

    public CNADrivers(DriverGenePanel panel) {
        this.oncoGenes = panel.oncoGenes();
        this.tsGenes = panel.tsGenes();
        this.deletionTargets = panel.deletionTargets();
        this.amplificationTargets = panel.amplificationTargets();
    }

    @NotNull
    public List<DriverCatalog> amplifications(final double ploidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        return geneCopyNumbers.stream()
                .filter(x -> x.minCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE)
                .filter(x -> oncoGenes.contains(x.gene()) | amplificationTargets.contains(x.gene()))
                .map(x -> ImmutableDriverCatalog.builder()
                        .chromosome(x.chromosome())
                        .chromosomeBand(x.chromosomeBand())
                        .gene(x.gene())
                        .missense(0)
                        .nonsense(0)
                        .inframe(0)
                        .frameshift(0)
                        .splice(0)
                        .dndsLikelihood(0)
                        .driverLikelihood(1)
                        .driver(DriverType.AMP)
                        .likelihoodMethod(LikelihoodMethod.AMP)
                        .category(tsGenes.contains(x.gene()) ? DriverCategory.TSG : DriverCategory.ONCO)
                        .biallelic(false)
                        .minCopyNumber(x.minCopyNumber())
                        .maxCopyNumber(x.maxCopyNumber())
                        .build())
                .collect(Collectors.toList());
    }

    @NotNull
    public List<DriverCatalog> deletions(@NotNull final List<GeneCopyNumber> geneCopyNumbers) {

        return geneCopyNumbers.stream()
                .filter(x -> x.minCopyNumber() < MAX_COPY_NUMBER_DEL)
                .filter(x -> tsGenes.contains(x.gene()) | deletionTargets.contains(x.gene()))
                .filter(x -> x.germlineHet2HomRegions() == 0 && x.germlineHomRegions() == 0)
                .map(x -> ImmutableDriverCatalog.builder()
                        .chromosome(x.chromosome())
                        .chromosomeBand(x.chromosomeBand())
                        .gene(x.gene())
                        .missense(0)
                        .nonsense(0)
                        .inframe(0)
                        .frameshift(0)
                        .splice(0)
                        .dndsLikelihood(0)
                        .driverLikelihood(1)
                        .driver(DriverType.DEL)
                        .likelihoodMethod(LikelihoodMethod.DEL)
                        .category(oncoGenes.contains(x.gene()) ? DriverCategory.ONCO : DriverCategory.TSG)
                        .biallelic(true)
                        .minCopyNumber(x.minCopyNumber())
                        .maxCopyNumber(x.maxCopyNumber())
                        .build())
                .collect(Collectors.toList());
    }
}
