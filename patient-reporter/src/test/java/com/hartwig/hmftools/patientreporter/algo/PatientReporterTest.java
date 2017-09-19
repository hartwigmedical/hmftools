package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusions;
import com.hartwig.hmftools.patientreporter.variants.StructuralVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.NullAnnotator;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String COSMIC_CSV = Resources.getResource("csv").getPath() + File.separator + "cosmic_gene_fusions.csv";

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.asMap());
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(geneModel);
        final StructuralVariantAnalyzer structuralVariantAnalyzer =
                new StructuralVariantAnalyzer(NullAnnotator.make(), geneModel.hmfRegions(), COSMICGeneFusions.readFromCSV(COSMIC_CSV));
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