package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.algo.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

final class SequencedReportDataLoader {

    private SequencedReportDataLoader() {
    }

    @NotNull
    static SequencedReportData buildFromFiles(@NotNull String fusionPairsLocation, @NotNull String promiscuousFiveLocation,
            @NotNull String promiscuousThreeLocation, @NotNull String drupGeneCsv, @NotNull String hotspotTsv,
            @NotNull String fastaFileLocation, @NotNull String highConfidenceBed) throws IOException {
        final GeneModel panelGeneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneList());
        final CompoundEnrichment compoundEnrichment = new CompoundEnrichment(HotspotEnrichment.fromHotspotsFile(hotspotTsv));

        final KnownFusionsModel knownFusionsModel = KnownFusionsModel.fromInputStreams(new FileInputStream(fusionPairsLocation),
                new FileInputStream(promiscuousFiveLocation),
                new FileInputStream(promiscuousThreeLocation));
        final DrupActionabilityModel drupActionabilityModel = new DrupActionabilityModel(drupGeneCsv);

        return ImmutableSequencedReportData.of(panelGeneModel,
                compoundEnrichment,
                knownFusionsModel,
                drupActionabilityModel,
                new IndexedFastaSequenceFile(new File(fastaFileLocation)),
                BEDFileLoader.fromBedFile(highConfidenceBed));
    }
}
