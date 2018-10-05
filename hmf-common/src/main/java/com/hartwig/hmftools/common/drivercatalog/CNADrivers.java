package com.hartwig.hmftools.common.drivercatalog;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

public class CNADrivers {

    private static final double MIN_COPY_NUMBER_RELATIVE_INCREASE = 3;
    private static final double MAX_COPY_NUMBER_DEL = 0.5;

    private final Set<String> oncoGenes;
    private final Set<String> tsGenes;
    private final Set<String> amplificationTargets;
    private final Set<String> deletionTargets;

    @NotNull
    public static Set<String> amplificationTargets() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/AmplificationTargets.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    @NotNull
    public static Set<String> deletionTargets() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/DeletionTargets.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    public CNADrivers() {
        this.amplificationTargets = amplificationTargets();
        this.oncoGenes = DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet();

        this.deletionTargets = deletionTargets();
        this.tsGenes = DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet();
    }

    @NotNull
    public List<DriverCatalog> amplifications(final double ploidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        return geneCopyNumbers.stream()
                .filter(x -> x.minCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE)
                .map(TranscriptRegion::gene)
                .filter(x -> oncoGenes.contains(x) | amplificationTargets.contains(x))
                .map(x -> ImmutableDriverCatalog.builder()
                        .gene(x)
                        .missense(0)
                        .nonsense(0)
                        .inframe(0)
                        .frameshift(0)
                        .splice(0)
                        .dndsLikelihood(0)
                        .driverLikelihood(1)
                        .driver(DriverType.AMP)
                        .category(tsGenes.contains(x) ? DriverCategory.TSG : DriverCategory.ONCO)
                        .build())
                .collect(Collectors.toList());
    }

    @NotNull
    public List<DriverCatalog> deletions(@NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        return geneCopyNumbers.stream()
                .filter(x -> x.minCopyNumber() < MAX_COPY_NUMBER_DEL)
                .filter(x -> x.germlineHet2HomRegions() == 0 && x.germlineHomRegions() == 0)
                .map(TranscriptRegion::gene)
                .filter(x -> tsGenes.contains(x) | deletionTargets.contains(x))
                .map(x -> ImmutableDriverCatalog.builder()
                        .gene(x)
                        .missense(0)
                        .nonsense(0)
                        .inframe(0)
                        .frameshift(0)
                        .splice(0)
                        .dndsLikelihood(0)
                        .driverLikelihood(1)
                        .driver(DriverType.DEL)
                        .category(oncoGenes.contains(x) ? DriverCategory.ONCO : DriverCategory.TSG)
                        .build())
                .collect(Collectors.toList());
    }
}
