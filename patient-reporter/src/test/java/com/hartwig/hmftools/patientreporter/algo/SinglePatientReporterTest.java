package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.slicing.HmfSlicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class SinglePatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String BED_DIRECTORY = Resources.getResource("bed").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final String hmfSlicingBed = BED_DIRECTORY + File.separator + "hmf_gene_panel.tsv";
        final HmfSlicer hmfSlicingRegion = SlicerFactory.fromHmfGenePanelFile(hmfSlicingBed);

        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion, hmfSlicingRegion, hmfSlicingRegion);
        final CopyNumberAnalyzer copyNumberAnalyzer = CopyNumberAnalyzer.fromHmfSlicingRegion(hmfSlicingRegion);

        final SinglePatientReporter algo =
                new SinglePatientReporter(false, buildTestCpctEcrfModel(), Lims.buildEmptyModel(), variantAnalyzer, copyNumberAnalyzer,
                        null);

        assertNotNull(algo.run(RUN_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
    }
}