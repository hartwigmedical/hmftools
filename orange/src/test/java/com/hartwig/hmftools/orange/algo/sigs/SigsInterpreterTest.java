package com.hartwig.hmftools.orange.algo.sigs;

import static org.junit.Assert.assertNotNull;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;

import org.junit.Test;

public class SigsInterpreterTest
{
    private static final List<SignatureAllocation> signature =
            Lists.newArrayList(ImmutableSignatureAllocation.builder().signature("Sig1").allocation(0D).percent(0D).build());
    private static final Map<String, String> etiologyPerSignature = new HashMap<>();

    static
    {
        etiologyPerSignature.put("Sig1", "APOBEC");
    }

    @Test
    public void shouldConvertAndAnnotateSignatures()
    {
        assertNotNull(SigsInterpreter.convertAndAnnotateWithEtiology(signature, etiologyPerSignature));
    }
}
