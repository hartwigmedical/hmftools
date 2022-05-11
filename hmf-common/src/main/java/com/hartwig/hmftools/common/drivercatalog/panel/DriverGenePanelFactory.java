package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.InputStream;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodFile;

import org.jetbrains.annotations.NotNull;

public final class DriverGenePanelFactory
{
    @NotNull
    public static DriverGenePanel empty()
    {
        return create(Collections.emptyList());
    }

    @NotNull
    public static DriverGenePanel create(final List<DriverGene> genes)
    {
        Map<String, DndsDriverGeneLikelihood> tsgLikelihood = tsgLikelihood().stream()
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        Map<String, DndsDriverGeneLikelihood> oncoLikelihood = oncoLikelihood().stream()
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        return ImmutableDriverGenePanel.builder()
                .driverGenes(genes)
                .tsgLikelihood(tsgLikelihood)
                .oncoLikelihood(oncoLikelihood)
                .build();
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> tsgLikelihood()
    {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/dnds_driver_likelihood_tsg.tsv"));
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> oncoLikelihood()
    {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/dnds_driver_likelihood_onco.tsv"));
    }

    @NotNull
    private static InputStream resource(@NotNull final String resource)
    {
        return DriverGenePanelFactory.class.getResourceAsStream(resource);
    }
}
