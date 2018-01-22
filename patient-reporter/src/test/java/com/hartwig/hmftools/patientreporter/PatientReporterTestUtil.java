package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGeneModel;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGenes;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.algo.MSIAnalyzer;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

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
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfGeneMap());
        final DrupFilter drupFilter = new DrupFilter(drupFilterPath);
        final CosmicGeneModel cosmicGeneModel = CosmicGenes.buildModelFromCsv(cosmicPath);
        final CosmicFusionModel fusionModel = CosmicFusions.readFromCSV(fusionPath);
        final MSIAnalyzer msiAnalyzer = testMSIAnalyzer();
        return ImmutableHmfReporterData.of(geneModel, cosmicGeneModel, drupFilter, fusionModel, msiAnalyzer);
    }

    @NotNull
    private static MSIAnalyzer testMSIAnalyzer() throws FileNotFoundException {
        return new MSIAnalyzer() {
            @Override
            @Nullable
            public IndexedFastaSequenceFile reference() {
                return null;
            }

            @Override
            public double analyzeVariants(@NotNull final List<SomaticVariant> variants) throws FileNotFoundException {
                return 0.91;
            }
        };
    }

    @NotNull
    public static BaseReporterData testBaseReporterData() throws IOException, EmptyFileException {
        final String centerPath = Resources.getResource("center").getPath() + File.separator + "centers.csv";
        final CpctEcrfModel ecrfModel = new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
        final Lims lims = LimsFactory.empty();
        final CenterModel centerModel = Center.readFromCSV(centerPath);

        return ImmutableBaseReporterData.of(ecrfModel, lims, centerModel, SIGNATURE_PATH);
    }
}
