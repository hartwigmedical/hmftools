package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
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

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(hotspot),
                Lists.newArrayList(),
                Lists.newArrayList());

        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene("reportable")
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .build();

        ReportableVariant variantNonMatch = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene("non-match")
                .chromosome(chromosome)
                .position(position + 1)
                .ref(ref)
                .alt(alt)
                .build();

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
        assertEquals(variantMatch.gene(), reportedEvidence.gene());

        assertEquals(reportedEvidence.protectSources().size(), 1);
        ProtectSource protectSourceReportedEvidence = findBySource(reportedEvidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, protectSourceReportedEvidence.evidenceType());

        ProtectEvidence unreportedEvidence = findByGene(evidences, "unreported");
        assertFalse(unreportedEvidence.reported());
        assertEquals(unreportedMatch.gene(), unreportedEvidence.gene());

        assertEquals(unreportedEvidence.protectSources().size(), 1);
        ProtectSource protectSourceUnreportedEvidence = findBySource(unreportedEvidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION, protectSourceUnreportedEvidence.evidenceType());
    }

    @Test
    public void canDetermineCodingEffectOther() {
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT,
                VariantEvidence.extractCodingEffectOther(
                        "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
    }

    @Test
    public void canDetermineEvent() {
        ReportableVariant variantProtein = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .canonicalHgvsProteinImpact("p.?")
                .build();

        ReportableVariant variantCoding = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .canonicalHgvsProteinImpact("c.100 A>T")
                .build();

        ReportableVariant variantEffect = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .canonicalHgvsProteinImpact("splice_variant")
                .build();

        assertEquals("p.Gly83fs",
                VariantEvidence.determineEvent(false,
                        variantProtein,
                        "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
        assertEquals("c.246_247delCG",
                VariantEvidence.determineEvent(false,
                        variantCoding,
                        "ENST00000579755|c.246_247delCG||frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
        assertEquals("frameshift_variant",
                VariantEvidence.determineEvent(false,
                        variantEffect,
                        "ENST00000579755|||frameshift_variant|NONSENSE_OR_FRAMESHIFT"));

        assertEquals("p.?",
                VariantEvidence.determineEvent(true,
                        variantProtein,
                        "ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
        assertEquals("c.100 A>T",
                VariantEvidence.determineEvent(true,
                        variantCoding,
                        "ENST00000579755|c.246_247delCG||frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
        assertEquals("splice_variant",
                VariantEvidence.determineEvent(true,
                        variantEffect,
                        "ENST00000579755|||frameshift_variant|NONSENSE_OR_FRAMESHIFT"));
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

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(),
                Lists.newArrayList(range),
                Lists.newArrayList());

        ReportableVariant variantMatch = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene(gene)
                .chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantOutsideRange = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene("outside range")
                .chromosome(chromosome)
                .position(start - 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongGene = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene("wrong gene")
                .chromosome(chromosome)
                .position(start + 1)
                .gene("other gene")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();
        ReportableVariant variantWrongMutationType = ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .gene("wrong mutation type")
                .chromosome(chromosome)
                .position(start + 1)
                .gene(gene)
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        List<ReportableVariant> reportable =
                Lists.newArrayList(variantMatch, variantOutsideRange, variantWrongGene, variantWrongMutationType);
        List<ProtectEvidence> evidences = variantEvidence.evidence(reportable, Lists.newArrayList(), Lists.newArrayList());

        assertEquals(1, evidences.size());
        ProtectEvidence evidence = findByGene(evidences, gene);
        assertTrue(evidence.reported());

        assertEquals(evidence.protectSources().size(), 1);
        ProtectSource protectSourceUnreportedEvidence = findBySource(evidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.EXON_MUTATION, protectSourceUnreportedEvidence.evidenceType());
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

        VariantEvidence variantEvidence = new VariantEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
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

        assertEquals(actEvidence.protectSources().size(), 1);
        ProtectSource protectSourceActEvidence = findBySource(actEvidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.ACTIVATION, protectSourceActEvidence.evidenceType());

        ProtectEvidence inactEvidence = findByGene(evidences, inactivatedGene);
        assertFalse(inactEvidence.reported());

        assertEquals(inactEvidence.protectSources().size(), 1);
        ProtectSource protectSourceInacEvidence = findBySource(inactEvidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.INACTIVATION, protectSourceInacEvidence.evidenceType());
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
                .from(VariantTestFactory.createTestReportableVariant())
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .gene(gene)
                .driverLikelihood(driverLikelihood)
                .build();
    }

    @NotNull
    private static ProtectSource findBySource(@NotNull Set<ProtectSource> sources, @NotNull Knowledgebase source) {
        for (ProtectSource protectSource : sources) {
            if (protectSource.source() == source) {
                return protectSource;
            }
        }

        throw new IllegalStateException("Could not find evidence with source: " + source);
    }
}