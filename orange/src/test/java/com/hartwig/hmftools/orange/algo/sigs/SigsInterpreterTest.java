package com.hartwig.hmftools.orange.algo.sigs;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SigsInterpreterTest
{
    private static final Map<String, String> ETIOLOGY_PER_SIGNATURE = new HashMap<>();

    static
    {
        ETIOLOGY_PER_SIGNATURE.put("Sig1", "APOBEC");
    }

    @Test
    public void shouldConvertAndAnnotateSignatures()
    {
        List<com.hartwig.hmftools.datamodel.sigs.SignatureAllocation> sig1 =
                SigsInterpreter.interpret(Lists.newArrayList(createSignature("Sig1")), ETIOLOGY_PER_SIGNATURE);
        List<com.hartwig.hmftools.datamodel.sigs.SignatureAllocation> misalloc =
                SigsInterpreter.interpret(Lists.newArrayList(createSignature("MISALLOC")), ETIOLOGY_PER_SIGNATURE);
        assertEquals(sig1.iterator().next().etiology(), "APOBEC");
        assertEquals(misalloc.iterator().next().etiology(), "-");
        assertNull(SigsInterpreter.interpret(null, ETIOLOGY_PER_SIGNATURE));
    }

    @NotNull
    private static SignatureAllocation createSignature(@NotNull String signature)
    {
        return ImmutableSignatureAllocation.builder().signature(signature).allocation(0D).percent(0D).build();
    }
}
