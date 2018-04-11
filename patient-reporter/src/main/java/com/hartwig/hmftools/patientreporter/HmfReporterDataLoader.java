package com.hartwig.hmftools.patientreporter;

import java.io.FileInputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.cosmic.genes.CosmicGeneModel;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGenes;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.variants.ImmutableMicrosatelliteAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.MicrosatelliteAnalyzer;

import org.jetbrains.annotations.NotNull;

final class HmfReporterDataLoader {
    private HmfReporterDataLoader() {
    }

    @NotNull
    static HmfReporterData buildFromFiles(@NotNull final String cosmicGeneFile, @NotNull final String fusionPairsLocation,
            @NotNull final String promiscuousFiveLocation, @NotNull final String promiscuousThreeLocation,
            @NotNull final String drupFilterFile, @NotNull final String fastaFileLocation) throws IOException {
        final GeneModel panelGeneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneMap());
        final CosmicGeneModel cosmicGeneModel = CosmicGenes.readFromCSV(cosmicGeneFile);
        final KnownFusionsModel knownFusionsModel = KnownFusionsModel.fromInputStreams(new FileInputStream(fusionPairsLocation),
                new FileInputStream(promiscuousFiveLocation),
                new FileInputStream(promiscuousThreeLocation));
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);
        final MicrosatelliteAnalyzer microsatelliteAnalyzer = ImmutableMicrosatelliteAnalyzer.of(fastaFileLocation);
        return ImmutableHmfReporterData.of(panelGeneModel, cosmicGeneModel, knownFusionsModel, drupFilter, microsatelliteAnalyzer);
    }
}
