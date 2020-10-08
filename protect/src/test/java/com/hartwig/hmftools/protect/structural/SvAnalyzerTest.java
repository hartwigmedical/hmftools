package com.hartwig.hmftools.protect.structural;

import static com.hartwig.hmftools.protect.structural.SvAnalysisDatamodelTestFactory.createTestDisruptionBuilder;
import static com.hartwig.hmftools.protect.structural.SvAnalysisDatamodelTestFactory.createTestFusionBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.ProtectTestFactory;

import org.junit.Test;

public class SvAnalyzerTest {

    private static final String DISRUPTED_GENE = "GENE";

    @Test
    public void canAnalyzeFusionsDisruptions() {
        List<LinxFusion> testFusions = Lists.newArrayList(createTestFusionBuilder().geneStart("X").geneEnd("Y").build());
        List<LinxBreakend> testDisruptions =
                Lists.newArrayList(createTestDisruptionBuilder().gene(DISRUPTED_GENE).junctionCopyNumber(1.0).undisruptedCopyNumber(1.0).build());

        SvAnalysis analysis = SvAnalyzer.run(testFusions, testDisruptions, ProtectTestFactory.loadTestActionabilityAnalyzer(), null);

        assertEquals(1, analysis.reportableFusions().size());
        assertEquals(1, analysis.reportableDisruptions().size());
        assertEquals(0, analysis.evidenceItems().size());
    }
}
