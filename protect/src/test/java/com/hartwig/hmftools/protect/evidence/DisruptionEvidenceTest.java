package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.KnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.protect.ServeTestFactory;
import com.hartwig.serve.datamodel.Knowledgebase;
import com.hartwig.serve.datamodel.gene.ActionableGene;
import com.hartwig.serve.datamodel.gene.GeneLevelEvent;
import com.hartwig.serve.datamodel.gene.ImmutableActionableGene;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DisruptionEvidenceTest {

    @Test
    public void canDetermineEvidenceForHomozygousDisruptions() {
        String geneAmp = "geneAmp";
        String geneInact = "geneInact";
        String geneDel = "geneDel";
        String geneUnder = "geneUnder";

        ActionableGene amp = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneAmp)
                .event(GeneLevelEvent.AMPLIFICATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene inactivation = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneInact)
                .event(GeneLevelEvent.INACTIVATION)
                .source(Knowledgebase.CKB)
                .build();
        ActionableGene deletion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneDel)
                .event(GeneLevelEvent.DELETION)
                .source(Knowledgebase.CKB)
                .build();

        ActionableGene underexpression = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneUnder)
                .event(GeneLevelEvent.UNDEREXPRESSION)
                .source(Knowledgebase.CKB)
                .build();

        DisruptionEvidence disruptionEvidence =
                new DisruptionEvidence(EvidenceTestFactory.create(), Lists.newArrayList(amp, inactivation, deletion, underexpression));

        HomozygousDisruption matchAmp = create(geneAmp);
        HomozygousDisruption matchInact = create(geneInact);
        HomozygousDisruption nonMatch = create("other gene");

        List<ProtectEvidence> evidences = disruptionEvidence.evidence(Lists.newArrayList(matchAmp, matchInact, nonMatch));

        assertEquals(1, evidences.size());
        ProtectEvidence evidence = evidences.get(0);
        assertTrue(evidence.reported());
        assertEquals(geneInact, evidence.gene());
        assertEquals(DisruptionEvidence.HOMOZYGOUS_DISRUPTION_EVENT, evidence.event());

        assertEquals(evidence.sources().size(), 1);
        assertEquals(EvidenceType.INACTIVATION, findByKnowledgebase(evidence.sources(), Knowledgebase.CKB).evidenceType());
    }

    @NotNull
    private static HomozygousDisruption create(@NotNull String gene) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    private static KnowledgebaseSource findByKnowledgebase(@NotNull Set<KnowledgebaseSource> sources, @NotNull Knowledgebase knowledgebaseToFind) {
        for (KnowledgebaseSource source : sources) {
            if (source.name() == knowledgebaseToFind) {
                return source;
            }
        }

        throw new IllegalStateException("Could not find evidence from source: " + knowledgebaseToFind);
    }
}