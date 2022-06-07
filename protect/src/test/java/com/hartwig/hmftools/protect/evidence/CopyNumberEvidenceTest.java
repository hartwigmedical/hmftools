package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberEvidenceTest {

    @Test
    public void canDetermineCopyNumberEvidence() {
        String geneAmp = "geneAmp";
        String geneDel = "geneDel";
        ActionableGene amp = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneAmp)
                .event(GeneLevelEvent.AMPLIFICATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene inactivation = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneDel)
                .event(GeneLevelEvent.INACTIVATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene fusion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneAmp)
                .event(GeneLevelEvent.FUSION)
                .source(Knowledgebase.CKB)
                .build();

        CopyNumberEvidence copyNumberEvidence =
                new CopyNumberEvidence(EvidenceTestFactory.create(), Lists.newArrayList(amp, inactivation, fusion));

        GainLoss reportableAmp = GainLossTestFactory.createGainLoss(geneAmp, CopyNumberInterpretation.FULL_GAIN);
        GainLoss reportableDel = GainLossTestFactory.createGainLoss(geneDel, CopyNumberInterpretation.FULL_LOSS);
        GainLoss ampOnOtherGene = GainLossTestFactory.createGainLoss("other gene", CopyNumberInterpretation.PARTIAL_GAIN);

        List<GainLoss> reportableGainLosses = Lists.newArrayList(reportableAmp, reportableDel, ampOnOtherGene);
        List<GainLoss> unreportedGainLosses = Lists.newArrayList();
        List<ProtectEvidence> evidences = copyNumberEvidence.evidence(reportableGainLosses, unreportedGainLosses);

        assertEquals(2, evidences.size());

        ProtectEvidence ampEvidence = find(evidences, geneAmp);
        assertTrue(ampEvidence.reported());
        assertEquals(ampEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.AMPLIFICATION, findByKnowledgebase(ampEvidence.sources(), Knowledgebase.CKB).evidenceType());

        ProtectEvidence delEvidence = find(evidences, geneDel);
        assertTrue(delEvidence.reported());
        assertEquals(delEvidence.sources().size(), 1);
        assertEquals(ProtectEvidenceType.INACTIVATION, findByKnowledgebase(delEvidence.sources(), Knowledgebase.CKB).evidenceType());
    }

    @NotNull
    private static ProtectEvidence find(@NotNull List<ProtectEvidence> evidences, @NotNull String geneToFind) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.gene().equals(geneToFind)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence for gene: " + geneToFind);
    }

    @NotNull
    private static ProtectSource findByKnowledgebase(@NotNull Set<ProtectSource> sources, @NotNull Knowledgebase knowledgebaseToFind) {
        for (ProtectSource source : sources) {
            if (source.name() == knowledgebaseToFind) {
                return source;
            }
        }

        throw new IllegalStateException("Could not find evidence from source: " + knowledgebaseToFind);
    }
}