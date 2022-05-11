package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.serve.Knowledgebase;
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
        int position = 10;
        String ref = "C";
        String alt = "T";

        ActionableHotspot hotspot = ImmutableActionableHotspot.builder()
                .from(ServeTestFactory.createTestActionableHotspot())
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .source(Knowledgebase.CKB)
                .build();

        VariantEvidence variantEvidence =
                new VariantEvidence(EvidenceTestFactory.create(), Lists.newArrayList(hotspot), Lists.newArrayList(), Lists.newArrayList());

        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene("reportable")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        ReportableVariant variantWrongPosition = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene("wrong position")
                .chromosome(chromosome)
                .position(position + 1)
                .ref(ref)
                .alt(alt)
                .build();

        ReportableVariant variantWrongAlt = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene("wrong alt")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt("G")
                .build();

        SomaticVariant unreportedMatch = SomaticVariantTestBuilderFactory.create()
                .gene("unreported")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        List<ProtectEvidence> evidences = variantEvidence.evidence(Lists.newArrayList(variantMatch, variantWrongAlt, variantWrongPosition),
                Lists.newArrayList(),
                Lists.newArrayList(unreportedMatch));

        assertEquals(2, evidences.size());

        ProtectEvidence reportedEvidence = findByGene(evidences, "reportable");
        assertTrue(reportedEvidence.reported());
        assertEquals(reportedEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, reportedEvidence.sources().iterator().next().evidenceType());

        ProtectEvidence unreportedEvidence = findByGene(evidences, "unreported");
        assertFalse(unreportedEvidence.reported());
        assertEquals(unreportedEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, unreportedEvidence.sources().iterator().next().evidenceType());
    }

    @Test
    public void canDetermineVariantEvidenceForRanges() {
        String gene = "gene";
        String chromosome = "1";
        int start = 5;
        int end = 15;
        MutationTypeFilter mutationTypeFilter = MutationTypeFilter.MISSENSE;

        ActionableRange range = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene(gene)
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .mutationType(mutationTypeFilter)
                .source(Knowledgebase.CKB)
                .build();

        VariantEvidence variantEvidence =
                new VariantEvidence(EvidenceTestFactory.create(), Lists.newArrayList(), Lists.newArrayList(range), Lists.newArrayList());

        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene(gene)
                .chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalHgvsCodingImpact("match")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantOutsideRange = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene(gene)
                .chromosome(chromosome)
                .position(start - 1)
                .gene(gene)
                .canonicalHgvsCodingImpact("outside range")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongGene = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene("wrong gene")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("wrong gene")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongMutationType = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .gene(gene)
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("wrong mutation type")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        List<ReportableVariant> reportable =
                Lists.newArrayList(variantMatch, variantOutsideRange, variantWrongGene, variantWrongMutationType);
        List<ProtectEvidence> evidences = variantEvidence.evidence(reportable, Lists.newArrayList(), Lists.newArrayList());

        assertEquals(1, evidences.size());

        ProtectEvidence evidence = findByGene(evidences, gene);
        assertTrue(evidence.reported());
        assertEquals("match", evidence.event());
        assertEquals(evidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.EXON_MUTATION, evidence.sources().iterator().next().evidenceType());
    }

    @Test
    public void canDetermineVariantEvidenceForGenes() {
        String activatedGene = "gene1";
        String inactivatedGene = "gene2";
        String amplifiedGene = "gene3";

        ActionableGene actionableGene1 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(activatedGene)
                .event(GeneLevelEvent.ACTIVATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene actionableGene2 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(inactivatedGene)
                .event(GeneLevelEvent.INACTIVATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene actionableGene3 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(amplifiedGene)
                .event(GeneLevelEvent.AMPLIFICATION)
                .source(Knowledgebase.CKB)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.create(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(actionableGene1, actionableGene2, actionableGene3));

        ReportableVariant driverOnActivatedGene = withGeneAndDriverLikelihood(activatedGene, 1D);
        ReportableVariant passengerOnInactivatedGene = withGeneAndDriverLikelihood(inactivatedGene, 0D);
        ReportableVariant driverOnAmplifiedGene = withGeneAndDriverLikelihood(amplifiedGene, 0D);
        ReportableVariant driverOnOtherGene = withGeneAndDriverLikelihood("other", 1D);

        List<ReportableVariant> reportableVariants =
                Lists.newArrayList(driverOnActivatedGene, passengerOnInactivatedGene, driverOnAmplifiedGene, driverOnOtherGene);

        List<SomaticVariant> unreportedVariants = Lists.newArrayList(SomaticVariantTestBuilderFactory.create()
                .reported(false)
                .gene(activatedGene)
                .canonicalCodingEffect(CodingEffect.NONE)
                .build());

        List<ProtectEvidence> evidences = variantEvidence.evidence(reportableVariants, Lists.newArrayList(), unreportedVariants);

        assertEquals(2, evidences.size());

        ProtectEvidence actEvidence = findByGene(evidences, activatedGene);
        assertTrue(actEvidence.reported());
        assertEquals(actEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.ACTIVATION, actEvidence.sources().iterator().next().evidenceType());

        ProtectEvidence inactEvidence = findByGene(evidences, inactivatedGene);
        assertFalse(inactEvidence.reported());
        assertEquals(inactEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.INACTIVATION, inactEvidence.sources().iterator().next().evidenceType());
    }

    @NotNull
    private static ProtectEvidence findByGene(@NotNull List<ProtectEvidence> evidences, @NotNull String geneToFind) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.gene().equals(geneToFind)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence for gene: " + geneToFind);
    }

    @NotNull
    private static ReportableVariant withGeneAndDriverLikelihood(@NotNull String gene, double driverLikelihood) {
        return ImmutableReportableVariant.builder()
                .from(VariantTestFactory.create())
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .gene(gene)
                .driverLikelihood(driverLikelihood)
                .build();
    }
}