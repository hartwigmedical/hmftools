package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cosmic.Cosmic;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.common.slicing.HmfSlicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.variants.StructuralVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.NullAnnotator;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String BED_DIRECTORY = Resources.getResource("bed").getPath();
    private static final String COSMIC_EXAMPLE_FILE = Resources.getResource("csv").getPath() + File.separator + "cosmic_slice.csv";

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final String hmfSlicingBed = BED_DIRECTORY + File.separator + "hmf_gene_panel.tsv";
        final HmfSlicer hmfSlicingRegion = SlicerFactory.fromHmfGenePanelFile(hmfSlicingBed);
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion, hmfSlicingRegion, hmfSlicingRegion);
        final StructuralVariantAnalyzer structuralVariantAnalyzer =
                new StructuralVariantAnalyzer(NullAnnotator.make(), hmfSlicingRegion, Cosmic.buildModelFromCsv(COSMIC_EXAMPLE_FILE));
        final PatientReporter algo =
                new PatientReporter(buildTestCpctEcrfModel(), LimsJsonModel.buildEmptyModel(), variantAnalyzer, structuralVariantAnalyzer);
        assertNotNull(algo.run(RUN_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
    }
}