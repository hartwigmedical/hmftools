package com.hartwig.hmftools.purple;

import java.io.InputStream;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihoodFile;
import com.hartwig.hmftools.common.driver.panel.DriverGene;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverGeneResource
{
    public final List<DriverGene> DriverGeneList;
    public final Map<String,DriverGene> DriverGeneMap;

    public final Map<String,DndsDriverGeneLikelihood> TsgLikelihoodMap;
    public final Map<String,DndsDriverGeneLikelihood> OncoLikelihoodMap;

    public DriverGeneResource()
    {
        DriverGeneList = Collections.emptyList();
        DriverGeneMap = Collections.emptyMap();
        TsgLikelihoodMap = Collections.emptyMap();
        OncoLikelihoodMap = Collections.emptyMap();
    }

    public DriverGeneResource(final List<DriverGene> driveGenes)
    {
        DriverGeneList = driveGenes;

        DriverGeneMap = Maps.newHashMap();
        for(DriverGene driverGene : driveGenes)
        {
            DriverGeneMap.put(driverGene.gene(), driverGene);
        }

        List<DndsDriverGeneLikelihood> tsgLikelihoods = DndsDriverGeneLikelihoodFile.fromInputStream(
                resource("/dnds/dnds_driver_likelihood_tsg.tsv"));

        TsgLikelihoodMap = tsgLikelihoods.stream().collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        List<DndsDriverGeneLikelihood> oncoLikelihoods = DndsDriverGeneLikelihoodFile.fromInputStream(
                resource("/dnds/dnds_driver_likelihood_onco.tsv"));

        OncoLikelihoodMap = oncoLikelihoods.stream().collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));
    }

    private static InputStream resource(@NotNull final String resource)
    {
        return DriverGeneResource.class.getResourceAsStream(resource);
    }

    /*
    public Set<DriverGene> amplificationTargets()
    {
        return targets(DriverGene::reportAmplification);
    }

    public Set<DriverGene> deletionTargets()
    {
        return targets(DriverGene::reportDeletion);
    }

    private Set<DriverGene> targets(@NotNull Predicate<DriverGene> filter)
    {
        return driverGenes().stream().filter(filter).collect(Collectors.toSet());
    }
    */
}
