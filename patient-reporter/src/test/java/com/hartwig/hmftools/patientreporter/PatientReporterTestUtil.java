package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.algo.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.algo.DrupActionabilityModelFactory;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;
import com.hartwig.hmftools.patientreporter.algo.GeneModelFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class PatientReporterTestUtil {

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature.png").getPath();
    private static final String RVA_LOGO = Resources.getResource("rvalogo/rva_logo_test.jpg").getPath();

    private static final String REF_GENOME_PATH = Resources.getResource("refgenome/ref.fasta").getPath();

    private static final String CENTER_CSV = Resources.getResource("center/centers.csv").getPath();

    private static final String KNOWLEDGEBASE_PATH = Resources.getResource("actionability").getPath();
    private static final String DRUP_GENES_CSV = Resources.getResource("csv/drup_genes.csv").getPath();
    private static final String HOTSPOT_TSV = Resources.getResource("csv/hotspots.tsv").getPath();
    private static final String FUSION_PAIRS_CSV = Resources.getResource("csv/fusion_pairs.csv").getPath();
    private static final String PROMISCUOUS_FIVE_CSV = Resources.getResource("csv/promiscuous_five.csv").getPath();
    private static final String PROMISCUOUS_THREE_CSV = Resources.getResource("csv/promiscuous_three.csv").getPath();

    private PatientReporterTestUtil() {
    }

    @NotNull
    public static DrupActionabilityModel testDrupActionabilityModel() throws IOException {
        return DrupActionabilityModelFactory.buildFromCsv(DRUP_GENES_CSV);
    }

    @NotNull
    public static SequencedReportData testSequencedReportData() throws IOException {
        DrupActionabilityModel drupActionabilityModel = testDrupActionabilityModel();
        GeneModel geneModel = GeneModelFactory.create(drupActionabilityModel);
        CompoundEnrichment compoundEnrichment = new CompoundEnrichment(HotspotEnrichment.fromHotspotsFile(HOTSPOT_TSV));

        return ImmutableSequencedReportData.of(geneModel,
                ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_PATH),
                compoundEnrichment,
                testKnownFusionModel(),
                new IndexedFastaSequenceFile(new File(REF_GENOME_PATH)),
                TreeMultimap.create());
    }

    @NotNull
    public static KnownFusionsModel testKnownFusionModel() throws IOException {
        return KnownFusionsModel.fromInputStreams(new FileInputStream(FUSION_PAIRS_CSV),
                new FileInputStream(PROMISCUOUS_FIVE_CSV),
                new FileInputStream(PROMISCUOUS_THREE_CSV));
    }

    @NotNull
    public static BaseReportData testBaseReportData() throws IOException {
        final List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
        final Lims lims = LimsFactory.empty();
        final CenterModel centerModel = Center.readFromCSV(CENTER_CSV);
        return ImmutableBaseReportData.of(patientTumorLocations, lims, centerModel, SIGNATURE_PATH, RVA_LOGO);
    }
}
