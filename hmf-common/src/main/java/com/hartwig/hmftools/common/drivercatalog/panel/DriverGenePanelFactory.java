package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.InputStream;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodFile;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

public final class DriverGenePanelFactory
{
    private DriverGenePanelFactory()
    {
    }

    @NotNull
    public static DriverGenePanel empty()
    {
        return create(RefGenomeVersion.V37, Collections.emptyList());
    }

    @NotNull
    public static DriverGenePanel create(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final List<DriverGene> genes)
    {
        final Map<String, String> dndsTsGenes = Maps.newHashMap();
        final Map<String, String> dndsOncoGenes = Maps.newHashMap();

        final List<DndsDriverGeneLikelihood> oncoLikelihoodList = oncoLikelihood();
        final Set<String> dndsGenes = oncoLikelihoodList.stream().map(DndsDriverGeneLikelihood::gene).collect(Collectors.toSet());
        final DndsGeneName dndsGeneName = new DndsGeneName(refGenomeVersion, dndsGenes);

        for(DriverGene gene : genes)
        {
            if(gene.reportMissenseAndInframe() || gene.reportNonsenseAndFrameshift() || gene.reportSplice())
            {
                boolean isValidGene = dndsGeneName.isValid(gene);

                if(!isValidGene)
                    continue;

                final String dndsName = dndsGeneName.dndsGeneName(gene);
                if(gene.likelihoodType().equals(DriverCategory.TSG))
                {
                    dndsTsGenes.put(dndsName, gene.gene());
                }
                else
                {
                    dndsOncoGenes.put(dndsName, gene.gene());
                }
            }
        }

        Map<String, DndsDriverGeneLikelihood> tsgLikelihood = tsgLikelihood().stream()
                .filter(x -> dndsTsGenes.containsKey(x.gene()))
                .map(x -> ImmutableDndsDriverGeneLikelihood.builder().from(x).gene(dndsTsGenes.get(x.gene())).build())
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        Map<String, DndsDriverGeneLikelihood> oncoLikelihood = oncoLikelihoodList.stream()
                .filter(x -> dndsOncoGenes.containsKey(x.gene()))
                .map(x -> ImmutableDndsDriverGeneLikelihood.builder().from(x).gene(dndsOncoGenes.get(x.gene())).build())
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, x -> x));

        return ImmutableDriverGenePanel.builder().driverGenes(genes).tsgLikelihood(tsgLikelihood).oncoLikelihood(oncoLikelihood).build();
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> tsgLikelihood()
    {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/DndsDriverLikelihoodTsg.tsv"));
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> oncoLikelihood()
    {
        return DndsDriverGeneLikelihoodFile.fromInputStream(resource("/dnds/DndsDriverLikelihoodOnco.tsv"));
    }

    @NotNull
    private static InputStream resource(@NotNull final String resource)
    {
        return DriverGenePanelFactory.class.getResourceAsStream(resource);
    }
}
