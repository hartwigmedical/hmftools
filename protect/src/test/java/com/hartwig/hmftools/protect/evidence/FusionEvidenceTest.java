package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.treatment.ImmutableTreatment;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionEvidenceTest {

    @Test
    public void canDetermineFusionEvidence() {
        ActionableGene promiscuous3_1 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("EGFR")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment1").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene promiscuous3_2 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("TP53")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment2").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene promiscuous3_3 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("PTEN")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment3").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene promiscuous5 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("BRAF")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment4").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene amp = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("KRAS")
                .event(GeneLevelEvent.AMPLIFICATION)
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment5").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusion())
                .geneUp("EML4")
                .geneDown("ALK")
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment6").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene activation = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("CDK4")
                .event(GeneLevelEvent.ACTIVATION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment7").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene any = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("DUX8")
                .event(GeneLevelEvent.ANY_MUTATION)
                .source(Knowledgebase.ACTIN)
                .treatment(ImmutableTreatment.builder().treament("treatment8").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene other = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("AB")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment9").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableFusion ig_pair = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusion())
                .geneUp("IGH")
                .geneDown("BCL2")
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment10").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene ig_fusion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("IGH")
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment11").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();
        ActionableGene ig_over = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("IGK")
                .event(GeneLevelEvent.OVEREXPRESSION)
                .source(Knowledgebase.CKB)
                .treatment(ImmutableTreatment.builder().treament("treatment12").drugClasses(Sets.newHashSet("drugClasses")).build())
                .build();

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(new KnownFusionData(KNOWN_PAIR, "EML4", "ALK", "", ""));
        knownFusionCache.addData(new KnownFusionData(PROMISCUOUS_5, "BRAF", "", "", ""));
        knownFusionCache.addData(new KnownFusionData(PROMISCUOUS_3, "", "EGFR", "", ""));
        knownFusionCache.addData(new KnownFusionData(PROMISCUOUS_3, "", "PTEN", "", ""));
        knownFusionCache.addData(new KnownFusionData(PROMISCUOUS_3, "", "TP53", "", ""));
        knownFusionCache.addData(new KnownFusionData(PROMISCUOUS_3, "", "CDK4", "", ""));
        knownFusionCache.addData(new KnownFusionData(IG_PROMISCUOUS, "IGH", "", "", ""));
        knownFusionCache.addData(new KnownFusionData(IG_PROMISCUOUS, "IGK", "", "", ""));
        knownFusionCache.addData(new KnownFusionData(IG_KNOWN_PAIR, "IGH", "BCL2", "", ""));

        FusionEvidence fusionEvidence = new FusionEvidence(EvidenceTestFactory.create(),
                Lists.newArrayList(promiscuous3_1,
                        promiscuous3_2,
                        promiscuous3_3,
                        promiscuous5,
                        amp,
                        activation,
                        any,
                        other,
                        ig_over,
                        ig_fusion),
                Lists.newArrayList(fusion, ig_pair),
                knownFusionCache);

        LinxFusion reportedFusionMatch = create("EML4", "ALK", true, KNOWN_PAIR.toString());
        LinxFusion reportedFusionUnMatch = create("NRG1", "NRG1", true, EXON_DEL_DUP.toString());
        LinxFusion reportedPromiscuousMatch5 = create("BRAF", "other gene", true, PROMISCUOUS_5.toString());
        LinxFusion reportedPromiscuousMatch3 = create("other gene", "EGFR", true, PROMISCUOUS_3.toString());
        LinxFusion reportedPromiscuousUnMatch3 = create("other gene", "KRAS", true, PROMISCUOUS_3.toString());
        LinxFusion reportedPromiscuousNonMatch = create("other gene", "TP53", true, PROMISCUOUS_3.toString());
        LinxFusion unreportedPromiscuousMatch = create("other gene", "PTEN", false, PROMISCUOUS_3.toString());
        LinxFusion reportedPromiscuousMatch = create("other gene", "CDK4", true, PROMISCUOUS_3.toString());
        LinxFusion reportedOtherMatch = create("other gene", "AB", false, NONE.toString());
        LinxFusion reportedIgPromiscuous = create("IGH", "other gene", false, IG_PROMISCUOUS.toString());
        LinxFusion reportedIgPromiscuousOver = create("IGK", "other gene", false, IG_PROMISCUOUS.toString());
        LinxFusion reportedIgKnown = create("IGH", "BCL2", false, IG_KNOWN_PAIR.toString());

        List<ProtectEvidence> evidences = fusionEvidence.evidence(Lists.newArrayList(reportedFusionMatch,
                        reportedFusionUnMatch,
                        reportedPromiscuousMatch5,
                        reportedPromiscuousMatch3,
                        reportedPromiscuousUnMatch3,
                        reportedPromiscuousNonMatch,
                        reportedPromiscuousMatch,
                        reportedIgKnown),
                Lists.newArrayList(unreportedPromiscuousMatch, reportedOtherMatch, reportedIgPromiscuous, reportedIgPromiscuousOver));

        assertEquals(10, evidences.size());

        ProtectEvidence evidence1 = findByFusion(evidences, reportedFusionMatch);
        assertTrue(evidence1.reported());
        assertEquals(evidence1.sources().size(), 1);
        assertEquals(ProtectEvidenceType.FUSION_PAIR, findByKnowledgebase(evidence1, Knowledgebase.CKB, "treatment6").evidenceType());

        ProtectEvidence evidence2 = findByFusion(evidences, reportedPromiscuousMatch5);
        assertTrue(evidence2.reported());
        assertEquals(evidence2.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence2, Knowledgebase.ACTIN, "treatment4").evidenceType());

        ProtectEvidence evidence3 = findByFusion(evidences, reportedPromiscuousMatch3);
        assertTrue(evidence3.reported());
        assertEquals(evidence3.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence3, Knowledgebase.ACTIN, "treatment1").evidenceType());

        ProtectEvidence evidence4 = findByFusion(evidences, reportedPromiscuousNonMatch);
        assertTrue(evidence4.reported());
        assertEquals(evidence4.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence4, Knowledgebase.ACTIN, "treatment2").evidenceType());

        ProtectEvidence evidence5 = findByFusion(evidences, unreportedPromiscuousMatch);
        assertFalse(evidence5.reported());
        assertEquals(evidence5.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence5, Knowledgebase.ACTIN, "treatment3").evidenceType());

        ProtectEvidence evidence6 = findByFusion(evidences, reportedPromiscuousMatch);
        assertTrue(evidence6.reported());
        assertEquals(evidence6.sources().size(), 1);
        assertEquals(ProtectEvidenceType.ACTIVATION, findByKnowledgebase(evidence6, Knowledgebase.ACTIN, "treatment7").evidenceType());

        ProtectEvidence evidence7 = findByFusion(evidences, reportedOtherMatch);
        assertFalse(evidence7.reported());
        assertEquals(evidence7.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence7, Knowledgebase.CKB, "treatment9").evidenceType());

        ProtectEvidence evidence8 = findByFusion(evidences, reportedIgKnown);
        assertFalse(evidence8.reported());
        assertEquals(evidence8.sources().size(), 1);
        assertEquals(ProtectEvidenceType.FUSION_PAIR, findByKnowledgebase(evidence8, Knowledgebase.CKB, "treatment10").evidenceType());

        ProtectEvidence evidence9 = findByFusion(evidences, reportedIgPromiscuous);
        assertFalse(evidence9.reported());
        assertEquals(evidence9.sources().size(), 1);
        assertEquals(ProtectEvidenceType.PROMISCUOUS_FUSION,
                findByKnowledgebase(evidence9, Knowledgebase.CKB, "treatment11").evidenceType());

        ProtectEvidence evidence10 = findByFusion(evidences, reportedIgPromiscuousOver);
        assertFalse(evidence10.reported());
        assertEquals(evidence10.sources().size(), 1);
        assertEquals(ProtectEvidenceType.OVER_EXPRESSION, findByKnowledgebase(evidence10, Knowledgebase.CKB, "treatment12").evidenceType());
    }

    @NotNull
    private static ProtectEvidence findByFusion(@NotNull List<ProtectEvidence> evidences, @NotNull LinxFusion fusion) {
        String event = ProtectEventGenerator.fusionEvent(fusion);
        for (ProtectEvidence evidence : evidences) {
            if (evidence.event().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Cannot find evidence with fusion event: " + event);
    }

    @Test
    public void canCorrectlyFilterOnExonRange() {
        int minExonUp = 5;
        int maxExonUp = 7;
        int minExonDown = 2;
        int maxExonDown = 4;

        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusion())
                .geneUp("EML4")
                .minExonUp(minExonUp)
                .maxExonUp(maxExonUp)
                .geneDown("ALK")
                .minExonDown(minExonDown)
                .maxExonDown(maxExonDown)
                .build();

        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.addData(new KnownFusionData(KNOWN_PAIR, "EML4", "ALK", "", ""));

        FusionEvidence fusionEvidence =
                new FusionEvidence(EvidenceTestFactory.create(), Lists.newArrayList(), Lists.newArrayList(fusion), knownFusionCache);

        ImmutableLinxFusion.Builder builder = linxFusionBuilder("EML4", "ALK", true, KNOWN_PAIR.toString());

        List<LinxFusion> onMinRange = Lists.newArrayList(builder.fusedExonUp(minExonUp).fusedExonDown(minExonDown).build());
        assertEquals(1, fusionEvidence.evidence(onMinRange, Lists.newArrayList()).size());

        List<LinxFusion> onMaxRange = Lists.newArrayList(builder.fusedExonUp(maxExonUp).fusedExonDown(maxExonDown).build());
        assertEquals(1, fusionEvidence.evidence(onMaxRange, Lists.newArrayList()).size());

        List<LinxFusion> upGeneExonTooLow = Lists.newArrayList(builder.fusedExonUp(minExonUp - 1).fusedExonDown(minExonDown).build());
        assertEquals(0, fusionEvidence.evidence(upGeneExonTooLow, Lists.newArrayList()).size());

        List<LinxFusion> upGeneExonTooHigh = Lists.newArrayList(builder.fusedExonUp(maxExonUp + 1).fusedExonDown(minExonDown).build());
        assertEquals(0, fusionEvidence.evidence(upGeneExonTooHigh, Lists.newArrayList()).size());

        List<LinxFusion> downGeneExonTooLow = Lists.newArrayList(builder.fusedExonUp(minExonUp).fusedExonDown(minExonDown - 1).build());
        assertEquals(0, fusionEvidence.evidence(downGeneExonTooLow, Lists.newArrayList()).size());

        List<LinxFusion> downGeneExonTooHigh = Lists.newArrayList(builder.fusedExonUp(maxExonUp).fusedExonDown(maxExonDown + 1).build());
        assertEquals(0, fusionEvidence.evidence(downGeneExonTooHigh, Lists.newArrayList()).size());
    }

    @NotNull
    private static LinxFusion create(@NotNull String geneStart, @NotNull String geneEnd, boolean reported, @NotNull String reportType) {
        return linxFusionBuilder(geneStart, geneEnd, reported, reportType).build();
    }

    @NotNull
    private static ImmutableLinxFusion.Builder linxFusionBuilder(@NotNull String geneStart, @NotNull String geneEnd, boolean reported,
            @NotNull String reportType) {
        return ImmutableLinxFusion.builder()
                .from(LinxTestFactory.createMinimalTestFusion())
                .geneStart(geneStart)
                .geneEnd(geneEnd)
                .reported(reported)
                .reportedType(reportType);
    }

    @NotNull
    private static ProtectSource findByKnowledgebase(@NotNull ProtectEvidence evidence, @NotNull Knowledgebase knowledgebaseToFind,
            @NotNull String treatment) {
        if (evidence.treatment().equals(treatment)) {
            for (ProtectSource source : evidence.sources()) {
                if (source.name() == knowledgebaseToFind) {
                    return source;
                }
            }
        }

        throw new IllegalStateException("Could not find evidence from source: " + knowledgebaseToFind);
    }
}