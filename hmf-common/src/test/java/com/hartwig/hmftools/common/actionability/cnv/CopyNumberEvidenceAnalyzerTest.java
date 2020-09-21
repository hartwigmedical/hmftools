package com.hartwig.hmftools.common.actionability.cnv;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberEvidenceAnalyzerTest {

    @Test
    public void actionabilityWorksForCNVs() {
        ActionableCopyNumber actionabilityCNV = ImmutableActionableCopyNumber.builder()
                .gene("ERBB2")
                .type(CopyNumberType.AMPLIFICATION)
                .source("oncoKb")
                .reference("ERBB2:amp")
                .drug("Trastuzumab")
                .drugsType("EGFR inhibitor")
                .cancerType("Breast Cancer")
                .source("civic")
                .level("A")
                .response("Responsive")
                .build();

        CopyNumberEvidenceAnalyzer cnvAnalyzer = new CopyNumberEvidenceAnalyzer(Lists.newArrayList(actionabilityCNV));

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzerTestFactory.buildWithOneCancerTypeMapping("Skin Melanoma", "4159");

        assertEquals(1, cnvAnalyzer.evidenceForCopyNumber(createTestReportableGainsAndLosses("ERBB2"), 2D, "Breast", cancerTypeAnalyzer).size());
    }

    @NotNull
    private static ReportableGainLoss createTestReportableGainsAndLosses (@NotNull String gene){
        return createTestReportableGainLossBuilder().gene(gene).build();
    }

    @NotNull
    public static ImmutableReportableGainLoss.Builder createTestReportableGainLossBuilder() {
        return ImmutableReportableGainLoss.builder()
                .chromosome("1")
                .chromosomeBand("band")
                .gene(Strings.EMPTY)
                .copies(10)
                .interpretation(CopyNumberInterpretation.GAIN);
    }

}