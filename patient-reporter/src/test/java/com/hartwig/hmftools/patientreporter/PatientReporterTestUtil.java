package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.center.CenterModelFactory;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.actionability.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.actionability.DrupActionabilityModelFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModelFactory;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class PatientReporterTestUtil {

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();

    private static final String REF_GENOME_PATH = Resources.getResource("refgenome/ref.fasta").getPath();

    private static final String KNOWLEDGEBASE_PATH = Resources.getResource("actionability").getPath();

    private static final String CENTER_CSV = Resources.getResource("csv/centers.csv").getPath();
    private static final String CENTER_MANUAL_CSV = Resources.getResource("csv/manual_mapping.csv").getPath();

    private static final String DRUP_GENES_CSV = Resources.getResource("csv/drup_genes.csv").getPath();
    private static final String HOTSPOT_TSV = Resources.getResource("csv/hotspots.tsv").getPath();

    private static final String FUSION_FILE = Resources.getResource("test_run/svAnalysis/CPCT11111111T_fusions.csv").getPath();
    private static final String DISRUPTION_FILE = Resources.getResource("test_run/svAnalysis/CPCT11111111T_disruptions.csv").getPath();

    private PatientReporterTestUtil() {
    }

    @NotNull
    public static DrupActionabilityModel testDrupActionabilityModel() throws IOException {
        return DrupActionabilityModelFactory.buildFromCsv(DRUP_GENES_CSV);
    }

    @NotNull
    public static ActionabilityAnalyzer testActionabilityAnalyzer() throws IOException {
        return ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_PATH);
    }

    @NotNull
    static SvAnalyzer testSvAnalyzerModel() throws IOException {
        return SvAnalyzer.fromFiles(FUSION_FILE, DISRUPTION_FILE);
    }

    @NotNull
    public static BaseReportData testBaseReportData() {
        try {
            final List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
            final Lims lims = LimsFactory.empty();
            final CenterModel centerModel = CenterModelFactory.readFromCSV(CENTER_CSV, CENTER_MANUAL_CSV);
            return ImmutableBaseReportData.of(patientTumorLocations, lims, centerModel, SIGNATURE_PATH, RVA_LOGO_PATH);
        } catch (IOException exception) {
            throw new IllegalStateException("Could not generate test base reporter data: " + exception.getMessage());
        }
    }

    @NotNull
    public static SequencedReportData testSequencedReportData() {
        try {
            DrupActionabilityModel drupActionabilityModel = testDrupActionabilityModel();
            GeneModel geneModel = GeneModelFactory.create(drupActionabilityModel);
            CompoundEnrichment compoundEnrichment = new CompoundEnrichment(HotspotEnrichment.fromHotspotsFile(HOTSPOT_TSV));

            return ImmutableSequencedReportData.of(geneModel,
                    testActionabilityAnalyzer(),
                    compoundEnrichment,
                    new IndexedFastaSequenceFile(new File(REF_GENOME_PATH)),
                    TreeMultimap.create());
        } catch (IOException exception) {
            throw new IllegalStateException("Could not generate test sequenced report data: " + exception.getMessage());
        }
    }
}
