package com.hartwig.hmftools.orange.algo.protect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.junit.Test;

public class EvidenceSelectorTest {

    @Test
    public void canInterpretProtectData() {
        ProtectSource trialSource = ProtectTestFactory.createSource(EvidenceSelector.TRIAL_SOURCES.iterator().next());
        ProtectSource evidenceSource = ProtectTestFactory.createSource(Knowledgebase.VICC_CGI);

        assert !EvidenceSelector.TRIAL_SOURCES.contains(evidenceSource.name());

        ProtectEvidence reportableEvidence = ProtectTestFactory.builder().reported(true).addSources(evidenceSource).build();
        ProtectEvidence reportableTrial = ProtectTestFactory.builder().reported(true).addSources(trialSource).build();
        ProtectEvidence unreportedEvidence = ProtectTestFactory.builder().reported(false).addSources(evidenceSource).build();
        ProtectEvidence unreportedTrial = ProtectTestFactory.builder().reported(false).addSources(trialSource).build();

        List<ProtectEvidence> evidences = Lists.newArrayList(reportableEvidence, reportableTrial, unreportedEvidence, unreportedTrial);
        ProtectInterpretedData protect = ProtectInterpreter.interpret(evidences);

        assertEquals(1, protect.reportableEvidences().size());
        assertTrue(protect.reportableEvidences().contains(reportableEvidence));

        assertEquals(1, protect.reportableTrials().size());
        assertTrue(protect.reportableTrials().contains(reportableTrial));

        assertEquals(1, protect.unreportedEvidences().size());
        assertTrue(protect.unreportedEvidences().contains(unreportedEvidence));

        assertEquals(1, protect.unreportedTrials().size());
        assertTrue(protect.unreportedTrials().contains(unreportedTrial));
    }
}