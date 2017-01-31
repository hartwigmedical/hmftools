package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterAlgoTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String BED_DIRECTORY = Resources.getResource("bed").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final String bedFile = BED_DIRECTORY + File.separator + "HMF_Slicing.bed";
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        // KODU: Every region in the HMF slicing bed should be convertable to a slicing annotation!
        for (final GenomeRegion region : slicer.regions()) {
            assertNotNull(HMFSlicingAnnotation.fromGenomeRegion(region));
        }

        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(slicer, slicer, slicer);
        final CopyNumberAnalyzer copyNumberAnalyzer = CopyNumberAnalyzer.fromHmfSlicingRegion(slicer);
        new PatientReporterAlgo(RUN_DIRECTORY, buildTestCpctEcrfModel(), variantAnalyzer, copyNumberAnalyzer,
                null, null, false).run();
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(Lists.newArrayList(), Lists.newArrayList());
    }
}