package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testActionabilityAnalyzer;
import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestDisruptionBuilder;
import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestFusionBuilder;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SvAnalyzerTest {

    private static final String DISRUPTED_GENE = "TP53";

    @Test
    public void canAnalyzeFusionsDisruptions() throws IOException {
        SvAnalyzer testAnalyzer = buildTestAnalyzer();

        List<GeneCopyNumber> geneCopyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().minCopyNumber(2D).maxCopyNumber(2D).gene(DISRUPTED_GENE).build());
        SvAnalysis analysis = testAnalyzer.run(geneCopyNumbers, testActionabilityAnalyzer(), null);

        assertEquals(1, analysis.reportableFusions().size());
        assertEquals(1, analysis.reportableDisruptions().size());
        assertEquals(0, analysis.evidenceItems().size());
    }

    @NotNull
    private static SvAnalyzer buildTestAnalyzer() {
        List<ReportableGeneFusion> testFusions = Lists.newArrayList(
                createTestFusionBuilder().geneStart("X").geneEnd("Y").build());

        List<ReportableDisruption> testDisruptions = Lists.newArrayList(
                createTestDisruptionBuilder().gene(DISRUPTED_GENE).build());

        return new SvAnalyzer(testFusions, testDisruptions);
    }
}
