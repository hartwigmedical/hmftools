package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionEvidenceTest {

    @Test
    public void canDetermineFusionEvidence() {
        String geneUp = "geneUp";
        String geneDown = "geneDown";
        String genePromiscuous = "genePromiscuous";
        ActionableGene promiscuous = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(genePromiscuous)
                .event(GeneLevelEvent.FUSION)
                .build();
        ActionableGene amp = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(genePromiscuous)
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();
        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusion())
                .geneUp(geneUp)
                .geneDown(geneDown)
                .build();

        FusionEvidence fusionEvidence = new FusionEvidence(ProtectTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(promiscuous, amp),
                Lists.newArrayList(fusion));

        LinxFusion fusionMatch = linxFusionBuilder().geneStart(geneUp).geneEnd(geneDown).build();
        LinxFusion promiscuousMatch = linxFusionBuilder().geneStart(genePromiscuous).geneEnd("other gene").build();
        LinxFusion promiscuousNonMatch = linxFusionBuilder().geneStart("other gene").geneEnd(geneDown).build();

        List<ProtectEvidence> evidenceItems =
                fusionEvidence.evidence(Lists.newArrayList(fusionMatch, promiscuousMatch, promiscuousNonMatch));

        assertEquals(2, evidenceItems.size());

        assertTrue(evidenceItems.get(0).reported());
        assertEquals(fusionMatch.genomicEvent(), evidenceItems.get(0).genomicEvent());

        assertTrue(evidenceItems.get(1).reported());
        assertEquals(promiscuousMatch.genomicEvent(), evidenceItems.get(1).genomicEvent());
    }

    @Test
    public void canCorrectlyFilterOnExonRange() {
        String geneUp = "geneUp";
        String geneDown = "geneDown";
        int minExonUp = 5;
        int maxExonUp = 7;
        int minExonDown = 2;
        int maxExonDown = 4;

        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .from(ServeTestFactory.createTestActionableFusion())
                .geneUp(geneUp)
                .minExonUp(minExonUp)
                .maxExonUp(maxExonUp)
                .geneDown(geneDown)
                .minExonDown(minExonDown)
                .maxExonDown(maxExonDown)
                .build();

        FusionEvidence fusionEvidence =
                new FusionEvidence(ProtectTestFactory.createTestEvidenceFactory(), Lists.newArrayList(), Lists.newArrayList(fusion));

        ImmutableLinxFusion.Builder builder = linxFusionBuilder().geneStart(geneUp).geneEnd(geneDown);

        // On min range
        assertEquals(1,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(minExonUp).fusedExonDown(minExonDown).build())).size());

        // On max range
        assertEquals(1,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(maxExonUp).fusedExonDown(maxExonDown).build())).size());

        // Up gene exon too low
        assertEquals(0,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(minExonUp - 1).fusedExonDown(minExonDown).build())).size());

        // Up gene exon too high
        assertEquals(0,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(maxExonUp + 1).fusedExonDown(minExonDown).build())).size());

        // Down gene exon too low
        assertEquals(0,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(minExonUp).fusedExonDown(minExonDown - 1).build())).size());

        // Down gene exon too high
        assertEquals(0,
                fusionEvidence.evidence(Lists.newArrayList(builder.fusedExonUp(maxExonUp).fusedExonDown(maxExonDown + 1).build())).size());
    }

    @NotNull
    private static ImmutableLinxFusion.Builder linxFusionBuilder() {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(0)
                .threePrimeBreakendId(0)
                .name(Strings.EMPTY)
                .reported(true)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.INFRAME)
                .likelihood(FusionLikelihoodType.HIGH)
                .chainLength(0)
                .chainLinks(0)
                .chainTerminated(true)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(0)
                .skippedExonsDown(0)
                .fusedExonUp(0)
                .fusedExonDown(0)
                .geneStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(0D);
    }
}