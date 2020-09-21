package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestUtil {

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();
    private static final String COMPANY_LOGO_PATH = Resources.getResource("company_logo/hartwig_logo_test.jpg").getPath();

    private static final String KNOWLEDGEBASE_DIRECTORY = Resources.getResource("actionability").getPath();

    private static final String GERMLINE_GENES_REPORTING_CSV = Resources.getResource("germline/germline_genes_reporting.csv").getPath();
    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();

    private PatientReporterTestUtil() {
    }

    @NotNull
    public static ActionabilityAnalyzer testActionabilityAnalyzer() throws IOException {
        return ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY);
    }

    @NotNull
    public static ReportData testReportData() {
        List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
        Lims lims = LimsFactory.empty();

        return ImmutableQCFailReportData.of(patientTumorLocations, lims, SIGNATURE_PATH, RVA_LOGO_PATH, COMPANY_LOGO_PATH);
    }

    @NotNull
    public static AnalysedReportData testAnalysedReportData() {
        try {
            GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(GERMLINE_GENES_REPORTING_CSV);
            SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);

            return ImmutableAnalysedReportData.builder()
                    .from(testReportData())
                    .actionabilityAnalyzer(testActionabilityAnalyzer())
                    .germlineReportingModel(germlineReportingModel)
                    .summaryModel(summaryModel)
                    .build();
        } catch (IOException exception) {
            throw new IllegalStateException("Could not generate test analysed report data: " + exception.getMessage());
        }
    }
}
