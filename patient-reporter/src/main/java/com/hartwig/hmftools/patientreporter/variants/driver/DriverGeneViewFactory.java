package com.hartwig.hmftools.patientreporter.variants.driver;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodSupplier;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneViewFactory {

    private DriverGeneViewFactory() {
    }

    @NotNull
    public static DriverGeneView create() {
        return ImmutableDriverGeneView.builder()
                .oncoDriverGenes(DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet())
                .tsgDriverGenes(DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet())
                .build();
    }
}
