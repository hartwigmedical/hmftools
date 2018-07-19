package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.cosmic.CosmicGeneModel;
import com.hartwig.hmftools.common.cosmic.CosmicGenes;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.variants.ImmutableMicrosatelliteAnalyzer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

final class HmfReporterDataLoader {
    private HmfReporterDataLoader() {
    }

    @NotNull
    static HmfReporterData buildFromFiles(@NotNull String cosmicGeneFile, @NotNull String fusionPairsLocation,
            @NotNull String promiscuousFiveLocation, @NotNull String promiscuousThreeLocation, @NotNull String drupFilterFile,
            @NotNull String fastaFileLocation, @NotNull String highConfidenceBed) throws IOException {
        final GeneModel panelGeneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneList());
        final CosmicGeneModel cosmicGeneModel = CosmicGenes.readFromCSV(cosmicGeneFile);
        final KnownFusionsModel knownFusionsModel = KnownFusionsModel.fromInputStreams(new FileInputStream(fusionPairsLocation),
                new FileInputStream(promiscuousFiveLocation),
                new FileInputStream(promiscuousThreeLocation));
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);

        return ImmutableHmfReporterData.of(panelGeneModel,
                cosmicGeneModel,
                knownFusionsModel,
                drupFilter,
                ImmutableMicrosatelliteAnalyzer.of(new IndexedFastaSequenceFile(new File(fastaFileLocation))),
                BEDFileLoader.fromBedFile(highConfidenceBed));
    }
}
