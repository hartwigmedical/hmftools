package com.hartwig.hmftools.common.sigs;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class SignatureTestFactory
{
    @NotNull
    public static ImmutableSignatureAllocation.Builder builder()
    {
        return ImmutableSignatureAllocation.builder().signature("Sig1").allocation(0D).percent(0D);
    }

    @NotNull
    public static Map<String, String> createEtiologiesMap()
    {
        Map<String, String> etiologiesMap = Maps.newHashMap();
        etiologiesMap.put("Sig1", "APOBEC");
        return etiologiesMap;
    }
}
