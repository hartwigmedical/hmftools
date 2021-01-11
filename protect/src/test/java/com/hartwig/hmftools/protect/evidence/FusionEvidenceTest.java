package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.protect.evidence.ProtectEvidenceTestFactory.createTestBaseEvent;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
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
        ActionableGene promiscuous =
                ImmutableActionableGene.builder().from(createTestBaseEvent()).gene(geneUp).event(GeneLevelEvent.FUSION).build();
        ActionableGene amp =
                ImmutableActionableGene.builder().from(createTestBaseEvent()).gene(geneDown).event(GeneLevelEvent.AMPLIFICATION).build();
        ActionableFusion fusion = ImmutableActionableFusion.builder().from(createTestBaseEvent()).geneUp(geneUp).geneDown(geneDown).build();

        FusionEvidence fusionEvidence = new FusionEvidence(Lists.newArrayList(promiscuous, amp), Lists.newArrayList(fusion));

        LinxFusion fusionMatch = linxFusionBuilder().geneStart(geneUp).geneEnd(geneDown).build();
        LinxFusion promiscuousMatch = linxFusionBuilder().geneStart(geneUp).geneEnd("other gene").build();
        LinxFusion promiscuousNonMatch = linxFusionBuilder().geneStart("other gene").geneEnd(geneDown).build();

        List<ProtectEvidence> evidenceItems =
                fusionEvidence.evidence(Sets.newHashSet(), Lists.newArrayList(fusionMatch, promiscuousMatch, promiscuousNonMatch));

        assertEquals(2, evidenceItems.size());

        assertTrue(evidenceItems.get(0).reported());
        assertEquals(fusionMatch.genomicEvent(), evidenceItems.get(0).genomicEvent());

        assertTrue(evidenceItems.get(1).reported());
        assertEquals(promiscuousMatch.genomicEvent(), evidenceItems.get(1).genomicEvent());
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