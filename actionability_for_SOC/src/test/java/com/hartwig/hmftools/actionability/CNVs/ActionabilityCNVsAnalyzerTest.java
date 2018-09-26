package com.hartwig.hmftools.actionability.CNVs;

import static org.junit.Assert.*;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.internal.$processor$.meta.$GsonMirrors;
import org.junit.Ignore;
import org.junit.Test;

public class ActionabilityCNVsAnalyzerTest {


    @Ignore
    public void ActionabilityWorksCNVs() {
        ActionabilityCNVs actionanilityCNV = ImmutableActionabilityCNVs.builder()
                .gene("ERBB2")
                .cnvType("Amplification")
                .source("oncoKb")
                .reference("ERBB2:amp")
                .drugsName("Trastuzumab")
                .drugsType("EGFR inhibitor")
                .cancerType("Breast Cancer")
                .levelSource("1")
                .hmfLevel("A")
                .evidenceType("Predictive")
                .significanceSource("Responsive")
                .hmfResponse("Responsive")
                .build();

        ActionabilityCNVsAnalyzer cnvAnalyzer = new ActionabilityCNVsAnalyzer(Lists.newArrayList(actionanilityCNV));

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

        assertTrue(cnvAnalyzer.actionableCNVs(geneCopyNumber, "Breast"));
        assertFalse(cnvAnalyzer.actionableCNVs(geneCopyNumber, "Skin"));
    }
}