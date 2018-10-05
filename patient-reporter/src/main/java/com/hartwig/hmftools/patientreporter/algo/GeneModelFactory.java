package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public final class GeneModelFactory {

    private GeneModelFactory() {
    }

    @NotNull
    public static GeneModel create(@NotNull DrupActionabilityModel drupActionabilityModel) {
        final List<HmfTranscriptRegion> somaticVariantGenePanel = Lists.newArrayList();
        final List<HmfTranscriptRegion> cnvGenePanel = Lists.newArrayList();
        final Map<String, DriverCategory> geneDriverCategoryMap = Maps.newHashMap();

        Map<String, HmfTranscriptRegion> geneMap = HmfGenePanelSupplier.allGenesMap();

        for (String oncoGene : DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet()) {
            HmfTranscriptRegion geneDefinition = fetchGeneDefinition(oncoGene, geneMap);
            somaticVariantGenePanel.add(geneDefinition);
            cnvGenePanel.add(geneDefinition);
            geneDriverCategoryMap.put(oncoGene, DriverCategory.ONCO);
        }

        for (String tsgGene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet()) {
            HmfTranscriptRegion geneDefinition = fetchGeneDefinition(tsgGene, geneMap);
            somaticVariantGenePanel.add(geneDefinition);
            cnvGenePanel.add(geneDefinition);
            geneDriverCategoryMap.put(tsgGene, DriverCategory.TSG);
        }

        for (String ampTarget : CNADrivers.amplificationTargets()) {
            cnvGenePanel.add(fetchGeneDefinition(ampTarget, geneMap));
            geneDriverCategoryMap.put(ampTarget, DriverCategory.ONCO);
        }

        for (String delTarget : CNADrivers.deletionTargets()) {
            cnvGenePanel.add(fetchGeneDefinition(delTarget, geneMap));
            geneDriverCategoryMap.put(delTarget, DriverCategory.TSG);
        }

        final Set<String> actionableDrupGenes = drupActionabilityModel.actionableGenes();
        for (String drupGene : actionableDrupGenes) {
            HmfTranscriptRegion geneDefinition = fetchGeneDefinition(drupGene, geneMap);
            somaticVariantGenePanel.add(geneDefinition);
            cnvGenePanel.add(geneDefinition);
        }

        return ImmutableGeneModel.of(somaticVariantGenePanel, cnvGenePanel, geneDriverCategoryMap, actionableDrupGenes);
    }

    @NotNull
    private static HmfTranscriptRegion fetchGeneDefinition(@NotNull String gene, @NotNull Map<String, HmfTranscriptRegion> geneMap) {
        HmfTranscriptRegion geneDefinition = geneMap.get(gene);
        if (geneDefinition == null) {
            throw new IllegalStateException("Every gene configured for the report should be defined in the gene map!");
        }
        return geneDefinition;
    }
}
