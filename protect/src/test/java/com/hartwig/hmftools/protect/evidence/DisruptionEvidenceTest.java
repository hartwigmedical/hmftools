package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.protect.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.linx.ReportableHomozygousDisruption;
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
        String gene = "gene";
        ActionableGene amp = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene)
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();
        ActionableGene inactivation = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene)
                .event(GeneLevelEvent.INACTIVATION)
                .build();
        ActionableGene deletion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene)
                .event(GeneLevelEvent.DELETION)
                .build();

        DisruptionEvidence disruptionEvidence =
                new DisruptionEvidence(ProtectTestFactory.createTestEvidenceFactory(), Lists.newArrayList(amp, inactivation, deletion));

        ReportableHomozygousDisruption match = create(gene);
        ReportableHomozygousDisruption nonMatch = create("other gene");

        List<ProtectEvidence> evidenceItems = disruptionEvidence.evidence(Lists.newArrayList(match, nonMatch));

        assertEquals(2, evidenceItems.size());

        assertTrue(evidenceItems.get(0).reported());
        assertTrue(evidenceItems.get(1).reported());
        assertEquals(match.genomicEvent(), evidenceItems.get(0).genomicEvent());
        assertEquals(match.genomicEvent(), evidenceItems.get(1).genomicEvent());
    }

    @NotNull
    private static ReportableHomozygousDisruption create(@NotNull String gene) {
        return ImmutableReportableHomozygousDisruption.builder().chromosome(Strings.EMPTY).chromosomeBand(Strings.EMPTY).gene(gene).build();
    }
}