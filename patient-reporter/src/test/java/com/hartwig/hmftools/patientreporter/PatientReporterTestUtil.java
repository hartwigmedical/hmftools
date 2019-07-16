package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModel;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModelFactory;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.hospital.HospitalModelFactory;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneViewFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class PatientReporterTestUtil {

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();
    private static final String COMPANY_LOGO_PATH = Resources.getResource("company_logo/hartwig_logo_test.jpg").getPath();

    private static final String KNOWLEDGEBASE_DIRECTORY = Resources.getResource("actionability").getPath();
    private static final String REF_GENOME_PATH = Resources.getResource("refgenome/ref.fasta").getPath();

    private static final String DRUP_GENES_CSV = Resources.getResource("actionability/drup_genes.csv").getPath();
    private static final String GERMLINE_GENES_REPORTING_CSV = Resources.getResource("csv/germline_genes_reporting.csv").getPath();
    private static final String SAMPLE_SUMMARY_CSV = Resources.getResource("csv/sample_summary.csv").getPath();

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
        HospitalModel hospitalModel = HospitalModelFactory.empty();

        return ImmutableQCFailReportData.of(patientTumorLocations, lims, hospitalModel, SIGNATURE_PATH, RVA_LOGO_PATH, COMPANY_LOGO_PATH);
    }

    @NotNull
    public static AnalysedReportData testAnalysedReportData() {
        try {
            DrupActionabilityModel drupActionabilityModel = DrupActionabilityModelFactory.buildFromCsv(DRUP_GENES_CSV);
            GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(GERMLINE_GENES_REPORTING_CSV);
            SummaryModel summaryModel = SummaryFile.buildFromCsv(SAMPLE_SUMMARY_CSV);

            return ImmutableAnalysedReportData.builder()
                    .from(testReportData())
                    .driverGeneView(DriverGeneViewFactory.create())
                    .drupActionabilityModel(drupActionabilityModel)
                    .actionabilityAnalyzer(testActionabilityAnalyzer())
                    .refGenomeFastaFile(new IndexedFastaSequenceFile(new File(REF_GENOME_PATH)))
                    .highConfidenceRegions(TreeMultimap.create())
                    .germlineReportingModel(germlineReportingModel)
                    .summaryModel(summaryModel)
                    .build();
        } catch (IOException exception) {
            throw new IllegalStateException("Could not generate test sequenced report data: " + exception.getMessage());
        }
    }
}
