package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.actionability.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.actionability.DrupActionabilityModelFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModelFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

final class SequencedReportDataLoader {

    private SequencedReportDataLoader() {
    }

    @NotNull
    static SequencedReportData buildFromFiles(@NotNull String knowledgebasePath, @NotNull String fusionPairsLocation,
            @NotNull String promiscuousFiveLocation, @NotNull String promiscuousThreeLocation, @NotNull String drupGeneCsv,
            @NotNull String hotspotTsv, @NotNull String fastaFileLocation, @NotNull String highConfidenceBed) throws IOException {
        final ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebasePath);

        final DrupActionabilityModel drupActionabilityModel = DrupActionabilityModelFactory.buildFromCsv(drupGeneCsv);
        final GeneModel panelGeneModel = GeneModelFactory.create(drupActionabilityModel);

        final KnownFusionsModel knownFusionsModel = KnownFusionsModel.fromInputStreams(new FileInputStream(fusionPairsLocation),
                new FileInputStream(promiscuousFiveLocation),
                new FileInputStream(promiscuousThreeLocation));

        return ImmutableSequencedReportData.of(panelGeneModel,
                actionabilityAnalyzer,
                HotspotEnrichment.fromHotspotsFile(hotspotTsv),
                knownFusionsModel,
                new IndexedFastaSequenceFile(new File(fastaFileLocation)),
                BEDFileLoader.fromBedFile(highConfidenceBed));
    }
}
