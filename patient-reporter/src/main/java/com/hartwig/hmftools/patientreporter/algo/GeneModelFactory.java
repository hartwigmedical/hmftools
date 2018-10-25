package com.hartwig.hmftools.patientreporter.algo;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.patientreporter.actionability.DrupActionabilityModel;

import org.jetbrains.annotations.NotNull;

public final class GeneModelFactory {

    private GeneModelFactory() {
    }

    @NotNull
    public static GeneModel create(@NotNull DrupActionabilityModel drupActionabilityModel) {
        Map<String, HmfTranscriptRegion> somaticVariantDriverGenePanel = Maps.newHashMap();
        Map<String, DriverCategory> geneDriverCategoryMap = Maps.newHashMap();

        Map<String, HmfTranscriptRegion> geneMap = HmfGenePanelSupplier.allGenesMap();

        for (String oncoGene : DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet()) {
            somaticVariantDriverGenePanel.put(oncoGene, fetchGeneDefinition(oncoGene, geneMap));
            geneDriverCategoryMap.put(oncoGene, DriverCategory.ONCO);
        }

        for (String tsgGene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet()) {
            somaticVariantDriverGenePanel.put(tsgGene, fetchGeneDefinition(tsgGene, geneMap));
            geneDriverCategoryMap.put(tsgGene, DriverCategory.TSG);
        }

        // KODU: For any gene that is actionable in DRUP and configured to onco/tsg, use this in driver category map,
        // only if not configured in the driver catalog.
        for (Map.Entry<String, DriverCategory> drupDriverCategoryEntry : drupActionabilityModel.geneDriverCategoryMap().entrySet()) {
            geneDriverCategoryMap.putIfAbsent(drupDriverCategoryEntry.getKey(), drupDriverCategoryEntry.getValue());
        }

        return ImmutableGeneModel.of(somaticVariantDriverGenePanel,
                toGeneMap(CNADrivers.amplificationTargets(), geneMap),
                toGeneMap(CNADrivers.deletionTargets(), geneMap),
                toGeneMap(drupActionabilityModel.actionableGenes(), geneMap),
                geneDriverCategoryMap);
    }

    @NotNull
    private static Map<String, HmfTranscriptRegion> toGeneMap(@NotNull Set<String> genes,
            @NotNull Map<String, HmfTranscriptRegion> geneMap) {
        Map<String, HmfTranscriptRegion> filteredGeneMap = Maps.newHashMap();
        for (String gene : genes) {
            filteredGeneMap.put(gene, fetchGeneDefinition(gene, geneMap));
        }
        return filteredGeneMap;
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
