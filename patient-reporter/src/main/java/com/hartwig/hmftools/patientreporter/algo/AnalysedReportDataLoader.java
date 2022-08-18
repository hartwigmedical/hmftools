package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.patientreporter.ReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.remarks.SpecialRemarkFile;
import com.hartwig.hmftools.patientreporter.remarks.SpecialRemarkModel;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;

import org.jetbrains.annotations.NotNull;

public final class AnalysedReportDataLoader {

    private AnalysedReportDataLoader() {
    }

    @NotNull
    public static AnalysedReportData buildFromFiles(@NotNull ReportData reportData, @NotNull String germlineReportingTsv,
             @NotNull String sampleSpecialRemarkTsv, @NotNull String knownFusionFile) throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(germlineReportingTsv);
        SpecialRemarkModel specialRemarkModel = SpecialRemarkFile.buildFromTsv(sampleSpecialRemarkTsv);

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        if (!knownFusionCache.loadFile(knownFusionFile)) {
            throw new IOException("Could not load known fusions from " + knownFusionFile);
        }

        return ImmutableAnalysedReportData.builder()
                .from(reportData)
                .germlineReportingModel(germlineReportingModel)
                .specialRemarkModel(specialRemarkModel)
                .knownFusionCache(knownFusionCache)
                .build();
    }
}
