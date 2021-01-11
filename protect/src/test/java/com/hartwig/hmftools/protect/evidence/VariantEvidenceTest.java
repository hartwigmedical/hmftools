package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantEvidenceTest {

    @Test
    public void canDetermineVariantEvidenceForHotspots() {
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

        ReportableVariant variantMatch =
                createTestReportableVariantBuilder().chromosome(chromosome).position(position).ref(ref).alt(alt).build();
        ReportableVariant variantNonMatch =
                createTestReportableVariantBuilder().chromosome(chromosome).position(position + 1).ref(ref).alt(alt).build();

        List<ProtectEvidenceItem> evidenceItems =
                variantEvidence.evidence(Sets.newHashSet(), Lists.newArrayList(variantMatch, variantNonMatch), Lists.newArrayList());

        assertEquals(1, evidenceItems.size());
        assertTrue(evidenceItems.get(0).reported());
        assertEquals(variantMatch.genomicEvent(), evidenceItems.get(0).genomicEvent());
    }

    @Test
    public void canDetermineVariantEvidenceForRanges() {
        String chromosome = "1";
        long start = 5;
        long end = 15;
        String gene = "gene";
        MutationTypeFilter mutationTypeFilter = MutationTypeFilter.MISSENSE;

        ActionableRange range = ImmutableActionableRange.builder()
                .from(createTestBaseEvent())
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .gene(gene)
                .mutationType(mutationTypeFilter)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(Lists.newArrayList(), Lists.newArrayList(range), Lists.newArrayList());

        ReportableVariant variantMatch = createTestReportableVariantBuilder().chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantOutsideRange = createTestReportableVariantBuilder().chromosome(chromosome)
                .position(start - 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongGene = createTestReportableVariantBuilder().chromosome(chromosome)
                .position(start + 1)
                .gene("other gene")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongMutationType = createTestReportableVariantBuilder().chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        List<ProtectEvidenceItem> evidenceItems = variantEvidence.evidence(Sets.newHashSet(),
                Lists.newArrayList(variantMatch, variantOutsideRange, variantWrongGene, variantWrongMutationType),
                Lists.newArrayList());

        assertEquals(1, evidenceItems.size());
        assertTrue(evidenceItems.get(0).reported());
        assertEquals(variantMatch.genomicEvent(), evidenceItems.get(0).genomicEvent());
    }

    @Test
    public void canDetermineVariantEvidenceForGenes() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        String gene3 = "gene3";

        ActionableGene actionableGene1 =
                ImmutableActionableGene.builder().from(createTestBaseEvent()).gene(gene1).event(GeneLevelEvent.ACTIVATION).build();
        ActionableGene actionableGene2 =
                ImmutableActionableGene.builder().from(createTestBaseEvent()).gene(gene2).event(GeneLevelEvent.INACTIVATION).build();
        ActionableGene actionableGene3 =
                ImmutableActionableGene.builder().from(createTestBaseEvent()).gene(gene3).event(GeneLevelEvent.AMPLIFICATION).build();

        VariantEvidence variantEvidence = new VariantEvidence(Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(actionableGene1, actionableGene2, actionableGene3));

        ReportableVariant variantMatchGene1 = createTestReportableVariantBuilder().gene(gene1).driverLikelihood(1D).build();
        ReportableVariant variantLowDriverGene2 = createTestReportableVariantBuilder().gene(gene2).driverLikelihood(0D).build();
        ReportableVariant variantMatchGene3 = createTestReportableVariantBuilder().gene(gene3).driverLikelihood(0D).build();
        ReportableVariant variantOtherGene = createTestReportableVariantBuilder().gene("other gene").driverLikelihood(1D).build();

        List<ProtectEvidenceItem> evidenceItems = variantEvidence.evidence(Sets.newHashSet(),
                Lists.newArrayList(variantMatchGene1, variantLowDriverGene2, variantMatchGene3, variantOtherGene),
                Lists.newArrayList());

        assertEquals(2, evidenceItems.size());
        assertTrue(evidenceItems.get(0).reported());
        assertEquals(variantMatchGene1.genomicEvent(), evidenceItems.get(0).genomicEvent());

        assertFalse(evidenceItems.get(1).reported());
        assertEquals(variantLowDriverGene2.genomicEvent(), evidenceItems.get(1).genomicEvent());
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