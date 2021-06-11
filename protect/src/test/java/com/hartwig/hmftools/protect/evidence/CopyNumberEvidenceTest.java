package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.protect.PurpleTestFactory;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.junit.Test;

public class CopyNumberEvidenceTest {

    @Test
    public void canDetermineCopyNumberEvidence() {
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
        ActionableGene fusion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(gene)
                .event(GeneLevelEvent.FUSION)
                .build();

        CopyNumberEvidence copyNumberEvidence =
                new CopyNumberEvidence(ProtectTestFactory.createTestEvidenceFactory(), Lists.newArrayList(amp, inactivation, fusion));

        ReportableGainLoss reportableAmp = PurpleTestFactory.testReportableGainLoss(gene, CopyNumberInterpretation.FULL_GAIN);
        ReportableGainLoss reportableDel = PurpleTestFactory.testReportableGainLoss(gene, CopyNumberInterpretation.FULL_LOSS);
        ReportableGainLoss ampOnOtherGene = PurpleTestFactory.testReportableGainLoss("other gene", CopyNumberInterpretation.PARTIAL_GAIN);

        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList(reportableAmp, reportableDel, ampOnOtherGene);
        List<ReportableGainLoss> unreportedGainLosses = Lists.newArrayList();
        List<ProtectEvidence> evidenceItems = copyNumberEvidence.evidence(reportableGainLosses, unreportedGainLosses);

        assertEquals(2, evidenceItems.size());

        assertTrue(evidenceItems.get(0).reported());
        assertEquals(reportableAmp.genomicEvent(), evidenceItems.get(0).genomicEvent());

        assertTrue(evidenceItems.get(1).reported());
        assertEquals(reportableDel.genomicEvent(), evidenceItems.get(1).genomicEvent());
    }
}