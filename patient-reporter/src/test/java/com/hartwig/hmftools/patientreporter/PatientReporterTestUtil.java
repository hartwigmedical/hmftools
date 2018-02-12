package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGeneModel;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGenes;
import com.hartwig.hmftools.common.ecrf.doid.TumorLocationDoidMapping;
import com.hartwig.hmftools.common.ecrf.projections.PatientCancerTypes;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.civic.AlterationAnalyzer;
import com.hartwig.hmftools.patientreporter.civic.CivicAnalyzer;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableAlteration;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableAlterationEvidence;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableAlterationMatch;
import com.hartwig.hmftools.patientreporter.variants.MicrosatelliteAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class PatientReporterTestUtil {

    public static final String SIGNATURE_PATH = Resources.getResource("signature").getPath() + File.separator + "signature.png";

    private PatientReporterTestUtil() {
    }

    @NotNull
    public static HmfReporterData testHmfReporterData() throws IOException, HartwigException {
        final String drupFilterPath = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
        final String cosmicPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_slice.csv";
        final String fusionPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_gene_fusions.csv";
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneMap());
        final CosmicGeneModel cosmicGeneModel = CosmicGenes.readFromCSV(cosmicPath);
        final CosmicFusionModel cosmicFusionModel = CosmicFusions.readFromCSV(fusionPath);
        final DrupFilter drupFilter = new DrupFilter(drupFilterPath);
        final MicrosatelliteAnalyzer microsatelliteAnalyzer = testMicrosatelliteAnalyzer();
        return ImmutableHmfReporterData.of(geneModel, cosmicGeneModel, cosmicFusionModel, drupFilter, microsatelliteAnalyzer);
    }

    @NotNull
    public static MicrosatelliteAnalyzer testMicrosatelliteAnalyzer() {
        return new MicrosatelliteAnalyzer() {
            @Override
            @Nullable
            public IndexedFastaSequenceFile reference() {
                return null;
            }

            @Override
            public double analyzeVariants(@NotNull final List<SomaticVariant> variants) {
                return 0.91;
            }
        };
    }

    @NotNull
    public static BaseReporterData testBaseReporterData() throws IOException, EmptyFileException {
        final String centerPath = Resources.getResource("center").getPath() + File.separator + "centers.csv";
        final List<PatientCancerTypes> patientsCancerTypes = Lists.newArrayList();
        final Lims lims = LimsFactory.empty();
        final CenterModel centerModel = Center.readFromCSV(centerPath);
        return ImmutableBaseReporterData.of(patientsCancerTypes, lims, centerModel, SIGNATURE_PATH);
    }

    @NotNull
    public static List<Alteration> mockedAlterations() {
        return Lists.newArrayList(ImmutableAlteration.of("BRAF",
                "p.Val600Glu",
                Lists.newArrayList(ImmutableAlterationEvidence.of("Sensitivity",
                        "Vemurafenib (B)\nDabrafenib (off-label: B)",
                        "CIViC",
                        "https://civic.genome.wustl.edu/events/genes/5/summary/variants/17/summary#variant")),
                Lists.newArrayList(ImmutableAlterationMatch.of("extended",
                        "V600",
                        "7",
                        "140453136",
                        "140453137",
                        "",
                        "",
                        "",
                        "protein_altering_variant",
                        "",
                        "https://civic.genome.wustl.edu/events/genes/5/summary/variants/17/summary#variant"))),
                ImmutableAlteration.of("ERBB2",
                        "GAIN",
                        Lists.newArrayList(ImmutableAlterationEvidence.of("Resistance or Non-Response",
                                "Erlotinib (C)\nGefitinib (C)\nPanitumumab (off-label:B)",
                                "CIViC",
                                "https://civic.genome.wustl.edu/events/genes/20/summary/variants/306/summary#variant")),
                        Lists.newArrayList(ImmutableAlterationMatch.of("manual",
                                "AMPLIFICATION",
                                "17",
                                "37856333",
                                "37884915",
                                "",
                                "",
                                "",
                                "transcript_amplification",
                                "",
                                "https://civic.genome.wustl.edu/events/genes/20/summary/variants/306/summary#variant"))));
    }

    @NotNull
    public static AlterationAnalyzer mockedCivicAnalyzer() {
        //@formatter:off
        return (@NotNull final List<VariantReport> reportedVariants, @NotNull final List<GeneCopyNumber> copyNumbers,
                @NotNull final List<GeneDisruptionData> disruptions, @NotNull final List<GeneFusionData> fusions, @NotNull final GeneModel geneModel,
                @NotNull final Set<String> tumorDoids) -> mockedAlterations();
        //@formatter:on
    }

    @SuppressWarnings("unused")
    @NotNull
    public static List<Alteration> runCivicAnalysis(List<VariantReport> variants, List<GeneCopyNumber> copyNumbers,
            @NotNull final List<GeneDisruptionData> disruptions, @NotNull final List<GeneFusionData> fusions, GeneModel geneModel,
            String tumorType) {
        final TumorLocationDoidMapping doidMapping = TumorLocationDoidMapping.fromResource("/tumor_location_doid_mapping.csv");
        return new CivicAnalyzer().run(variants, copyNumbers, disruptions, fusions, geneModel, doidMapping.doidsForTumorType(tumorType));
    }

}
