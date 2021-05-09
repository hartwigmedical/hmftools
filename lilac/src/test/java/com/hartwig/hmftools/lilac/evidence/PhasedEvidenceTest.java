package com.hartwig.hmftools.lilac.evidence;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.junit.Test;

public class PhasedEvidenceTest
{
    @Test
    public void testInconsistentEvidence()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("CAT", 4);
        evidence.put("ATC", 5);
        PhasedEvidence victim = new PhasedEvidence(Lists.newArrayList(new Integer(0), new Integer(1), new Integer(3)), evidence);
        HlaSequenceLoci catCandidate = create(new HlaSequence(HlaAllele.fromString("A*01:01"), "CART"));
        HlaSequenceLoci atcCandidate = create(new HlaSequence(HlaAllele.fromString("A*01:02"), "ATRC"));
        HlaSequenceLoci wildAtcCandidate = create(new HlaSequence(HlaAllele.fromString("A*01:03"), "*TRC"));
        HlaSequenceLoci wildCandidate = create(new HlaSequence(HlaAllele.fromString("A*01:04"), "****"));

        PhasedEvidence noMissing = victim.inconsistentEvidence(Lists.newArrayList(catCandidate, atcCandidate));
        assertEquals(0, noMissing.totalEvidence());

        PhasedEvidence catMissing = victim.inconsistentEvidence(Lists.newArrayList(atcCandidate));
        assertEquals(1, catMissing.getEvidence().size());
        assertTrue(catMissing.getEvidence().containsKey("CAT"));

        PhasedEvidence wildAtcMatch = victim.inconsistentEvidence(Lists.newArrayList(wildAtcCandidate));
        assertEquals(1, wildAtcMatch.getEvidence().size());
        assertTrue(wildAtcMatch.getEvidence().containsKey("CAT"));

        PhasedEvidence wildMatch = victim.inconsistentEvidence(Lists.newArrayList(wildCandidate));
        assertEquals(0, wildMatch.totalEvidence());
    }

    private HlaSequenceLoci create(final HlaSequence sequences)
    {
        return HlaSequenceLoci.create(sequences.Allele, sequences.getRawSequence(), sequences.getRawSequence());
    }

}
