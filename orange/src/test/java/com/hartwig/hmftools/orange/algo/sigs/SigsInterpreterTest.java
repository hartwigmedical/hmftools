package com.hartwig.hmftools.orange.algo.sigs;

import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_MISALLOCATED;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;

import org.junit.Test;

public class SigsInterpreterTest
{
    private static final Map<String, String> ETIOLOGY_PER_SIGNATURE = new HashMap<>();

    static
    {
        ETIOLOGY_PER_SIGNATURE.put("Sig1", "APOBEC");
    }

    @Test
    public void testSignatureConversion()
    {
        String sig1 = "Sig1";
        String sig2 = "Sig2";
        String sig3 = "Sig3";

        List<SignatureAllocation> signatures = List.of(
                ImmutableSignatureAllocation.builder().signature(sig1).allocation(7).percent(0.5).build(),
                ImmutableSignatureAllocation.builder().signature(sig2).allocation(20).percent(0.05).build(),
                ImmutableSignatureAllocation.builder().signature(sig3).allocation(8).percent(0.08).build(),
                ImmutableSignatureAllocation.builder().signature(SIG_MISALLOCATED).allocation(10).percent(0.1).build()
        );

        Map<String,String> etiologies = Map.of(sig1, sig1, sig2, sig2, sig3, sig3);
        List<com.hartwig.hmftools.datamodel.sigs.SignatureAllocation> convertedSigs = SigsInterpreter.interpret(signatures, etiologies);

        assertEquals(2, convertedSigs.size());
        assertEquals(sig3, convertedSigs.get(0).signature());
        assertEquals(SIG_MISALLOCATED, convertedSigs.get(1).signature());
    }
}
