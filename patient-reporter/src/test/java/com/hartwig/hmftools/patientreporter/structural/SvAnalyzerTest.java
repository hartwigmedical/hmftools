package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testActionabilityAnalyzer;
import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestDisruptionBuilder;
import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestFusionBuilder;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.ReportableDisruption;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;

import org.junit.Test;

public class SvAnalyzerTest {

    private static final String DISRUPTED_GENE = "GENE";

    @Test
    public void canAnalyzeFusionsDisruptions() throws IOException {
        List<ReportableGeneFusion> testFusions = Lists.newArrayList(createTestFusionBuilder().geneStart("X").geneEnd("Y").build());
        List<ReportableDisruption> testDisruptions =
                Lists.newArrayList(createTestDisruptionBuilder().gene(DISRUPTED_GENE).ploidy(1.0).undisruptedCopyNumber(1.0).build());

        SvAnalysis analysis = SvAnalyzer.run(testFusions, testDisruptions, testActionabilityAnalyzer(), null);

        assertEquals(1, analysis.reportableFusions().size());
        assertEquals(1, analysis.reportableDisruptions().size());
        assertEquals(0, analysis.evidenceItems().size());
    }
}
