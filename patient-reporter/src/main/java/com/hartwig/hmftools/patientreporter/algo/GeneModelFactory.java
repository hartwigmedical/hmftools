package com.hartwig.hmftools.patientreporter.algo;

import java.util.Map;

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
        final Map<String, HmfTranscriptRegion> somaticVariantDriverGenePanel = Maps.newHashMap();
        final Map<String, DriverCategory> somaticVariantDriverCategoryMap = Maps.newHashMap();

        final Map<String, HmfTranscriptRegion> significantlyAmplifiedGenes = Maps.newHashMap();
        final Map<String, HmfTranscriptRegion> significantlyDeletedGenes = Maps.newHashMap();

        Map<String, HmfTranscriptRegion> geneMap = HmfGenePanelSupplier.allGenesMap();

        for (String oncoGene : DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet()) {
            somaticVariantDriverGenePanel.put(oncoGene, fetchGeneDefinition(oncoGene, geneMap));
            somaticVariantDriverCategoryMap.put(oncoGene, DriverCategory.ONCO);
        }

        for (String tsgGene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet()) {
            somaticVariantDriverGenePanel.put(tsgGene, fetchGeneDefinition(tsgGene, geneMap));
            somaticVariantDriverCategoryMap.put(tsgGene, DriverCategory.TSG);
        }

        for (String ampTarget : CNADrivers.amplificationTargets()) {
            significantlyAmplifiedGenes.put(ampTarget, fetchGeneDefinition(ampTarget, geneMap));
        }

        for (String delTarget : CNADrivers.deletionTargets()) {
            significantlyDeletedGenes.put(delTarget, fetchGeneDefinition(delTarget, geneMap));
        }

        Map<String, HmfTranscriptRegion> drupActionabilityOPanel = Maps.newHashMap();
        for (String drupGene : drupActionabilityModel.actionableGenes()) {
            drupActionabilityOPanel.put(drupGene, fetchGeneDefinition(drupGene, geneMap));
        }

        return ImmutableGeneModel.of(somaticVariantDriverGenePanel,
                somaticVariantDriverCategoryMap,
                significantlyAmplifiedGenes,
                significantlyDeletedGenes,
                drupActionabilityOPanel);
    }

    @NotNull
    private static HmfTranscriptRegion fetchGeneDefinition(@NotNull String gene, @NotNull Map<String, HmfTranscriptRegion> geneMap) {
        HmfTranscriptRegion geneDefinition = geneMap.get(gene);
        if (geneDefinition == null) {
            throw new IllegalStateException("Every gene configured for the report should be defined in the gene map: " + gene);
        }
        return geneDefinition;
    }
}
