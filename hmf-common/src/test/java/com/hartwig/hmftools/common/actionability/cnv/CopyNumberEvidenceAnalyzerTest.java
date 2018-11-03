package com.hartwig.hmftools.common.actionability.cnv;

import static org.junit.Assert.assertFalse;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

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

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzerTestFactory.buildWithOneCancerTypeMapping("Skin Melanoma", 4159);

        GeneCopyNumber geneCopyNumber = ImmutableGeneCopyNumber.builder()
                .gene("ERBB2")
                .maxCopyNumber(45)
                .minCopyNumber(30.001)
                .somaticRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .minRegions(1)
                .minRegionStart(1)
                .minRegionEnd(10)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .nonsenseBiallelicCount(0)
                .nonsenseNonBiallelicCount(0)
                .nonsenseNonBiallelicPloidy(0)
                .spliceBiallelicCount(0)
                .spliceNonBiallelicCount(0)
                .missenseNonBiallelicPloidy(0)
                .minMinorAllelePloidy(0)
                .transcriptID("trans")
                .transcriptVersion(1)
                .chromosomeBand("12.1")
                .chromosome("1")
                .start(1)
                .end(2)
                .spliceNonBiallelicPloidy(0)
                .missenseBiallelicCount(0)
                .missenseNonBiallelicCount(0)
                .build();

           assertFalse(cnvAnalyzer.evidenceForCopyNumber(geneCopyNumber, "Breast", cancerTypeAnalyzer).isEmpty());
    }
}