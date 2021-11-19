package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantTestFactory;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

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
                .from(ServeTestFactory.createTestActionableHotspot())
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(hotspot),
                Lists.newArrayList(),
                Lists.newArrayList());

        ReportableVariant base = VariantTestFactory.createTestReportableVariant();
        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(base)
                .gene("reportable")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();
        ReportableVariant variantNonMatch =
                ImmutableReportableVariant.builder().from(base).chromosome(chromosome).position(position + 1).ref(ref).alt(alt).build();
        SomaticVariant unreportedMatch = SomaticVariantTestBuilderFactory.create()
                .gene("unreported")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        List<ProtectEvidence> evidences = variantEvidence.evidence(Lists.newArrayList(variantMatch, variantNonMatch),
                Lists.newArrayList(),
                Lists.newArrayList(unreportedMatch));

        assertEquals(2, evidences.size());
        ProtectEvidence reportedEvidence = findByGene(evidences, "reportable");
        assertTrue(reportedEvidence.reported());
        assertEquals(variantMatch.genomicEvent(), reportedEvidence.genomicEvent());

        ProtectEvidence unreportedEvidence = findByGene(evidences, "unreported");
        assertFalse(unreportedEvidence.reported());
        assertEquals(unreportedMatch.genomicEvent(), unreportedEvidence.genomicEvent());
    }

    @NotNull
    private static ProtectEvidence findByGene(@NotNull List<ProtectEvidence> evidences, @NotNull String gene) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.genomicEvent().contains(gene)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence based on gene: " + gene);
    }

    @Test
    public void canDetermineVariantEvidenceForRanges() {
        String chromosome = "1";
        long start = 5;
        long end = 15;
        String gene = "gene";
        MutationTypeFilter mutationTypeFilter = MutationTypeFilter.MISSENSE;

        ActionableRange range = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .gene(gene)
                .mutationType(mutationTypeFilter)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(),
                Lists.newArrayList(range),
                Lists.newArrayList());

        ReportableVariant base = VariantTestFactory.createTestReportableVariant();
        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(base)
                .chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantOutsideRange = ImmutableReportableVariant.builder()
                .from(base)
                .chromosome(chromosome)
                .position(start - 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongGene = ImmutableReportableVariant.builder()
                .from(base)
                .chromosome(chromosome)
                .position(start + 1)
                .gene("other gene")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongMutationType = ImmutableReportableVariant.builder()
                .from(base)
                .chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        List<ProtectEvidence> evidenceItems =
                variantEvidence.evidence(Lists.newArrayList(variantMatch, variantOutsideRange, variantWrongGene, variantWrongMutationType),
                        Lists.newArrayList(),
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

        ActionableGene actionableGene1 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene1)
                .event(GeneLevelEvent.ACTIVATION)
                .build();
        ActionableGene actionableGene2 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene2)
                .event(GeneLevelEvent.INACTIVATION)
                .build();
        ActionableGene actionableGene3 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene3)
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(actionableGene1, actionableGene2, actionableGene3));

        ReportableVariant base = VariantTestFactory.createTestReportableVariant();
        ReportableVariant variantMatchGene1 = ImmutableReportableVariant.builder().from(base).gene(gene1).driverLikelihood(1D).build();
        ReportableVariant variantLowDriverGene2 = ImmutableReportableVariant.builder().from(base).gene(gene2).driverLikelihood(0D).build();
        ReportableVariant variantMatchGene3 = ImmutableReportableVariant.builder().from(base).gene(gene3).driverLikelihood(0D).build();
        ReportableVariant variantOtherGene = ImmutableReportableVariant.builder().from(base).gene("other").driverLikelihood(1D).build();

        List<ProtectEvidence> evidenceItems =
                variantEvidence.evidence(Lists.newArrayList(variantMatchGene1, variantLowDriverGene2, variantMatchGene3, variantOtherGene),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(2, evidenceItems.size());
        assertTrue(evidenceItems.get(0).reported());
        assertEquals(variantMatchGene1.genomicEvent(), evidenceItems.get(0).genomicEvent());

        assertFalse(evidenceItems.get(1).reported());
        assertEquals(variantLowDriverGene2.genomicEvent(), evidenceItems.get(1).genomicEvent());
    }
}