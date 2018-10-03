package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.patientreporter.algo.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class PatientReporterTestUtil {

    private static final String SIGNATURE_PATH = Resources.getResource("signature").getPath() + File.separator + "signature.png";

    private static final String REF_GENOME_PATH = Resources.getResource("refgenome").getPath() + File.separator + "ref.fasta";

    private static final String DRUP_GENES_CSV = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
    private static final String HOTSPOT_TSV = Resources.getResource("csv").getPath() + File.separator + "hotspots.tsv";
    private static final String FUSION_PAIRS_CSV = Resources.getResource("csv").getPath() + File.separator + "fusion_pairs.csv";
    private static final String PROMISCUOUS_FIVE_CSV = Resources.getResource("csv").getPath() + File.separator + "promiscuous_five.csv";
    private static final String PROMISCUOUS_THREE_CSV = Resources.getResource("csv").getPath() + File.separator + "promiscuous_three.csv";

    private PatientReporterTestUtil() {
    }

    @NotNull
    public static SequencedReportData testSequencedReportData() throws IOException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneList());
        CompoundEnrichment compoundEnrichment = new CompoundEnrichment(HotspotEnrichment.fromHotspotsFile(HOTSPOT_TSV));
        final DrupActionabilityModel drupActionabilityModel = new DrupActionabilityModel(DRUP_GENES_CSV);

        return ImmutableSequencedReportData.of(geneModel,
                compoundEnrichment,
                testKnownFusionModel(),
                drupActionabilityModel,
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
        final String centerPath = Resources.getResource("center").getPath() + File.separator + "centers.csv";
        final List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
        final Lims lims = LimsFactory.empty();
        final CenterModel centerModel = Center.readFromCSV(centerPath);
        return ImmutableBaseReportData.of(patientTumorLocations, lims, centerModel, SIGNATURE_PATH);
    }
}
