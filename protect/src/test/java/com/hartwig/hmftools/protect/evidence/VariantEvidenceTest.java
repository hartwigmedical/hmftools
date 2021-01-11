package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.protect.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantEvidenceTest {

    @Test
    public void canDetermineVariantEvidence() {
        String chromosome = "1";
        long position = 10;
        String ref = "C";
        String alt = "T";

        ActionableHotspot hotspot = ImmutableActionableHotspot.builder()
                .from(createTestBaseEvent())
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(Lists.newArrayList(hotspot), Lists.newArrayList(), Lists.newArrayList());

        ReportableVariant variant =
                createTestReportableVariantBuilder().chromosome(chromosome).position(position).ref(ref).alt(alt).build();

        List<ProtectEvidenceItem> evidenceItems =  variantEvidence.evidence(Sets.newHashSet(), variant);

        assertEquals(1, evidenceItems.size());
    }

    @NotNull
    private static ActionableEvent createTestBaseEvent() {
        return ActionabilityTestUtil.create(Knowledgebase.VICC_CGI,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false);
    }
}