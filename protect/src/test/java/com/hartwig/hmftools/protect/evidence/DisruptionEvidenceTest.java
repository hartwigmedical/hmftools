package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DisruptionEvidenceTest {

    @Test
    public void canDetermineEvidenceForHomozygousDisruptions() {
        String geneAmp = "geneAmp";
        String geneInact = "geneInact";
        String geneDel = "geneDel";
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

        DisruptionEvidence disruptionEvidence =
                new DisruptionEvidence(EvidenceTestFactory.createTestEvidenceFactory(), Lists.newArrayList(amp, inactivation, deletion));

        ReportableHomozygousDisruption matchAmp = create(geneAmp);
        ReportableHomozygousDisruption matchInact = create(geneInact);
        ReportableHomozygousDisruption nonMatch = create("other gene");

        List<ProtectEvidence> evidences = disruptionEvidence.evidence(Lists.newArrayList(matchAmp, matchInact, nonMatch));

        assertEquals(1, evidences.size());
        ProtectEvidence evidence = evidences.get(0);
        assertTrue(evidence.reported());
        assertEquals(geneInact, evidence.gene());
        assertEquals(DisruptionEvidence.HOMOZYGOUS_DISRUPTION_EVENT, evidence.event());

        assertEquals(evidence.protectSources().size(), 1);
        ProtectSource protectSource = findBySource(evidence.protectSources(), Knowledgebase.CKB);
        assertEquals(ProtectEvidenceType.INACTIVATION, protectSource.evidenceType());
    }

    @NotNull
    private static ReportableHomozygousDisruption create(@NotNull String gene) {
        return ImmutableReportableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    private static ProtectSource findBySource(@NotNull Set<ProtectSource> sources, @NotNull Knowledgebase source) {
        for (ProtectSource protectSource : sources) {
            if (protectSource.sources() == source) {
                return protectSource;
            }
        }

        throw new IllegalStateException("Could not find evidence with source: " + source);
    }
}