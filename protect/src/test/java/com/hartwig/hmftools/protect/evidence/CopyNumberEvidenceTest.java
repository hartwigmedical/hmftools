package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
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
                .build();
        ActionableGene inactivation = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneDel)
                .event(GeneLevelEvent.INACTIVATION)
                .build();
        ActionableGene fusion = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene(geneAmp)
                .event(GeneLevelEvent.FUSION)
                .build();

        CopyNumberEvidence copyNumberEvidence =
                new CopyNumberEvidence(EvidenceTestFactory.createTestEvidenceFactory(), Lists.newArrayList(amp, inactivation, fusion));

        ReportableGainLoss reportableAmp = PurpleTestFactory.createReportableGainLoss(geneAmp, CopyNumberInterpretation.FULL_GAIN);
        ReportableGainLoss reportableDel = PurpleTestFactory.createReportableGainLoss(geneDel, CopyNumberInterpretation.FULL_LOSS);
        ReportableGainLoss ampOnOtherGene = PurpleTestFactory.createReportableGainLoss("other gene", CopyNumberInterpretation.PARTIAL_GAIN);

        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList(reportableAmp, reportableDel, ampOnOtherGene);
        List<ReportableGainLoss> unreportedGainLosses = Lists.newArrayList();
        List<ProtectEvidence> evidences = copyNumberEvidence.evidence(reportableGainLosses, unreportedGainLosses);

        assertEquals(2, evidences.size());

        ProtectEvidence ampEvidence = find(evidences, geneAmp);
        assertTrue(ampEvidence.reported());
        assertEquals(reportableAmp.gene(), ampEvidence.gene());
        assertEquals(ProtectEvidenceType.AMPLIFICATION, ampEvidence.evidenceType());

        ProtectEvidence delEvidence = find(evidences, geneDel);
        assertTrue(delEvidence.reported());
        assertEquals(reportableDel.gene(), delEvidence.gene());
        assertEquals(ProtectEvidenceType.INACTIVATION, delEvidence.evidenceType());
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
}