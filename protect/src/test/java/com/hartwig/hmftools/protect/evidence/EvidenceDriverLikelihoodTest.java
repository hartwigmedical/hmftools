package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.DriverInterpretation;

import org.junit.Test;

public class EvidenceDriverLikelihoodTest {

    @Test
    public void canInterpretVariant(){
        assertTrue(EvidenceDriverLikelihood.interpretVariant(DriverInterpretation.HIGH));
        assertFalse(EvidenceDriverLikelihood.interpretVariant(DriverInterpretation.MEDIUM));
        assertFalse(EvidenceDriverLikelihood.interpretVariant(DriverInterpretation.LOW));
    }

    @Test
    public void canInterpretSignatures() {
        assertNull(EvidenceDriverLikelihood.interpretSignatures());
    }

    @Test
    public void canInterpretFusions(){
        assertTrue(EvidenceDriverLikelihood.interpretFusions(FusionLikelihoodType.HIGH));
        assertFalse(EvidenceDriverLikelihood.interpretFusions(FusionLikelihoodType.LOW));
        assertFalse(EvidenceDriverLikelihood.interpretFusions(FusionLikelihoodType.NA));
    }

    @Test
    public void canInterpretDisruptions(){
        assertTrue(EvidenceDriverLikelihood.interpretDisruptions());
    }

    @Test
    public void canInterpretCNV(){
        assertTrue(EvidenceDriverLikelihood.interpretCNV());
    }

    @Test
    public void canInterpretChord(){
        assertNull(EvidenceDriverLikelihood.interpretChord());
    }

    @Test
    public void canInterpretVirus(){
        assertTrue(EvidenceDriverLikelihood.interpretVirus());
    }
}