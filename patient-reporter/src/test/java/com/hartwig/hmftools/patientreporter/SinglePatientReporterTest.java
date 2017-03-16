package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.GenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.lims.TumorPercentageTestFactory;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotationFactory;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class SinglePatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String BED_DIRECTORY = Resources.getResource("bed").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final String hmfSlicingBed = BED_DIRECTORY + File.separator + "HMF_Slicing.bed";
        final Slicer hmfSlicingRegion = SlicerFactory.fromBedFile(hmfSlicingBed);
        // KODU: Every region in the HMF slicing bed should be convertable to a slicing annotation!
        for (final GenomeRegion region : hmfSlicingRegion.regions()) {
            assertNotNull(HMFSlicingAnnotationFactory.fromGenomeRegion(region));
        }

        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion, hmfSlicingRegion,
                hmfSlicingRegion);
        final CopyNumberAnalyzer copyNumberAnalyzer = CopyNumberAnalyzer.fromHmfSlicingRegion(hmfSlicingRegion);

        final SinglePatientReporter algo = new SinglePatientReporter(buildTestCpctEcrfModel(), variantAnalyzer,
                copyNumberAnalyzer, TumorPercentageTestFactory.buildTestTumorPercentages(), null);

        assertNotNull(algo.run(RUN_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(Lists.newArrayList(), Lists.newArrayList());
    }
}