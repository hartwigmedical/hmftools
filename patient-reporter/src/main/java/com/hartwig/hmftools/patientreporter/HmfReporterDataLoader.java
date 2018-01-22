package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGeneModel;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGenes;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.algo.ImmutableMSIAnalyzer;
import com.hartwig.hmftools.patientreporter.algo.MSIAnalyzer;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.jetbrains.annotations.NotNull;

public final class HmfReporterDataLoader {
    private HmfReporterDataLoader() {
    }

    @NotNull
    public static HmfReporterData buildFromFiles(@NotNull final String drupFilterFile, @NotNull final String cosmicFile,
            @NotNull final String fusionFile, @NotNull final String fastaFileLocation) throws IOException, HartwigException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfGeneMap());
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);
        final CosmicGeneModel cosmicGeneModel = CosmicGenes.buildModelFromCsv(cosmicFile);
        final CosmicFusionModel fusionModel = CosmicFusions.readFromCSV(fusionFile);
        final MSIAnalyzer msiAnalyzer = ImmutableMSIAnalyzer.of(fastaFileLocation);
        return ImmutableHmfReporterData.of(geneModel, cosmicGeneModel, drupFilter, fusionModel, msiAnalyzer);
    }
}
