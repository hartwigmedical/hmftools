package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.ServeTestFactory;
import com.hartwig.hmftools.common.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.common.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.common.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.common.serve.datamodel.MutationTypeFilter;
import com.hartwig.hmftools.common.serve.datamodel.gene.GeneLevelEvent;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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
                .from(ReportableVariantTestFactory.create())
                .gene("reportable")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        ReportableVariant variantWrongPosition = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("wrong position")
                .chromosome(chromosome)
                .position(position + 1)
                .ref(ref)
                .alt(alt)
                .build();

        ReportableVariant variantWrongAlt = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("wrong alt")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt("G")
                .build();

        SomaticVariant unreportedMatch =
                SomaticVariantTestFactory.builder().gene("unreported").chromosome(chromosome).position(position).ref(ref).alt(alt).build();

        List<ProtectEvidence> evidences = variantEvidence.evidence(Lists.newArrayList(variantMatch, variantWrongAlt, variantWrongPosition),
                Lists.newArrayList(),
                Lists.newArrayList(unreportedMatch));

        assertEquals(2, evidences.size());

        ProtectEvidence reportedEvidence = findByGene(evidences, "reportable");
        assertTrue(reportedEvidence.reported());
        assertEquals(reportedEvidence.sources().size(), 1);
        assertEquals(EvidenceType.HOTSPOT_MUTATION, reportedEvidence.sources().iterator().next().evidenceType());

        ProtectEvidence unreportedEvidence = findByGene(evidences, "unreported");
        assertFalse(unreportedEvidence.reported());
        assertEquals(unreportedEvidence.sources().size(), 1);
        assertEquals(EvidenceType.HOTSPOT_MUTATION, unreportedEvidence.sources().iterator().next().evidenceType());
    }

    @Test
    public void canDetermineVariantEvidenceForRanges() {
        String chromosome = "1";
        int start = 5;
        int end = 15;
        MutationTypeFilter mutationTypeFilter = MutationTypeFilter.MISSENSE;

        ActionableRange rangeHigh = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene("geneHigh")
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .mutationType(mutationTypeFilter)
                .source(Knowledgebase.CKB)
                .build();

        ActionableRange rangeMedium = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene("geneMedium")
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .mutationType(mutationTypeFilter)
                .source(Knowledgebase.CKB)
                .build();

        ActionableRange rangeLow = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .gene("geneLow")
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .mutationType(mutationTypeFilter)
                .source(Knowledgebase.CKB)
                .build();

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.create(),
                Lists.newArrayList(),
                Lists.newArrayList(rangeHigh, rangeMedium, rangeLow),
                Lists.newArrayList());

        ReportableVariant variantMatchHigh = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("geneHigh")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("match")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .driverLikelihood(0.9)
                .build();
        ReportableVariant variantMatchMedium = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("geneMedium")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("match")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .driverLikelihood(0.5)
                .build();
        ReportableVariant variantMatchLow = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("geneLow")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("match")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .driverLikelihood(0.1)
                .build();
        ReportableVariant variantOutsideRange = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("gene")
                .chromosome(chromosome)
                .position(start - 1)
                .canonicalHgvsCodingImpact("outside range")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongGene = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("wrong gene")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("wrong gene")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongMutationType = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("gene")
                .chromosome(chromosome)
                .position(start + 1)
                .canonicalHgvsCodingImpact("wrong mutation type")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        List<ReportableVariant> reportable = Lists.newArrayList(variantMatchHigh,
                variantMatchMedium,
                variantMatchLow,
                variantOutsideRange,
                variantWrongGene,
                variantWrongMutationType);
        List<ProtectEvidence> evidences = variantEvidence.evidence(reportable, Lists.newArrayList(), Lists.newArrayList());

        assertEquals(3, evidences.size());

        ProtectEvidence evidenceHigh = findByGene(evidences, "geneHigh");
        assertTrue(evidenceHigh.reported());
        assertEquals("match", evidenceHigh.event());
        assertEquals(evidenceHigh.sources().size(), 1);
        assertEquals(EvidenceType.EXON_MUTATION, evidenceHigh.sources().iterator().next().evidenceType());

        ProtectEvidence evidenceMedium = findByGene(evidences, "geneMedium");
        assertFalse(evidenceMedium.reported());
        assertEquals("match", evidenceMedium.event());
        assertEquals(evidenceMedium.sources().size(), 1);
        assertEquals(EvidenceType.EXON_MUTATION, evidenceMedium.sources().iterator().next().evidenceType());

        ProtectEvidence evidenceLow = findByGene(evidences, "geneLow");
        assertFalse(evidenceLow.reported());
        assertEquals("match", evidenceLow.event());
        assertEquals(evidenceLow.sources().size(), 1);
        assertEquals(EvidenceType.EXON_MUTATION, evidenceLow.sources().iterator().next().evidenceType());
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

        List<SomaticVariant> unreportedVariants = Lists.newArrayList(SomaticVariantTestFactory.builder()
                .reported(false)
                .gene(activatedGene)
                .canonicalCodingEffect(CodingEffect.NONE)
                .build());

        List<ProtectEvidence> evidences = variantEvidence.evidence(reportableVariants, Lists.newArrayList(), unreportedVariants);

        assertEquals(2, evidences.size());

        ProtectEvidence actEvidence = findByGene(evidences, activatedGene);
        assertTrue(actEvidence.reported());
        assertEquals(actEvidence.sources().size(), 1);
        assertEquals(EvidenceType.ACTIVATION, actEvidence.sources().iterator().next().evidenceType());

        ProtectEvidence inactEvidence = findByGene(evidences, inactivatedGene);
        assertFalse(inactEvidence.reported());
        assertEquals(inactEvidence.sources().size(), 1);
        assertEquals(EvidenceType.INACTIVATION, inactEvidence.sources().iterator().next().evidenceType());
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
                .from(ReportableVariantTestFactory.create())
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .gene(gene)
                .driverLikelihood(driverLikelihood)
                .build();
    }
}