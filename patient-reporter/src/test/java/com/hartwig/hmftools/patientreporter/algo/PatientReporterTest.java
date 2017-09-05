package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.hmfslicer.HmfGeneRegionSupplier;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final GeneModel geneModel = new GeneModel(HmfGeneRegionSupplier.asMap());
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(geneModel, geneModel.slicer(), geneModel.slicer());
        final PatientReporter algo = new PatientReporter(buildTestCpctEcrfModel(), LimsJsonModel.buildEmptyModel(), variantAnalyzer);
        assertNotNull(algo.run(RUN_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
    }
}