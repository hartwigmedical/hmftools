package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModel;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModelFactory;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModelFactory;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

final class SequencedReportDataLoader {

    private SequencedReportDataLoader() {
    }

    @NotNull
    static SequencedReportData buildFromFiles(@NotNull String knowledgebaseDir, @NotNull String drupGeneCsv, @NotNull String hotspotTsv,
            @NotNull String fastaFileLocation, @NotNull String highConfidenceBed, @NotNull String germlineGenesCsv,
            @NotNull String summarySamplesCSV) throws IOException {
        final ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDir);

        final DrupActionabilityModel drupActionabilityModel = DrupActionabilityModelFactory.buildFromCsv(drupGeneCsv);
        final GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(germlineGenesCsv);
        final SummaryModel summaryModel = SummaryFile.buildFromCsv(summarySamplesCSV);
        final GeneModel panelGeneModel = GeneModelFactory.create(drupActionabilityModel);

        return ImmutableSequencedReportData.of(panelGeneModel,
                actionabilityAnalyzer,
                HotspotEnrichment.fromHotspotsFile(hotspotTsv),
                new IndexedFastaSequenceFile(new File(fastaFileLocation)),
                BEDFileLoader.fromBedFile(highConfidenceBed),
                germlineReportingModel,
                summaryModel);
    }
}
