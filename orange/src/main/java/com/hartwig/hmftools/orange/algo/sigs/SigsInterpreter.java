package com.hartwig.hmftools.orange.algo.sigs;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SigsInterpreter
{
    @Nullable
    public static List<SignatureAllocation> interpret(
            @Nullable List<com.hartwig.hmftools.common.sigs.SignatureAllocation> signatureAllocation,
            @NotNull Map<String, String> etiologyPerSignature)
    {
        if(signatureAllocation == null)
        {
            return null;
        }
        return signatureAllocation.stream()
                .map(signature -> ImmutableSignatureAllocation.builder()
                        .signature(signature.signature())
                        .etiology(annotateWithEtiology(signature, etiologyPerSignature))
                        .allocation(signature.allocation())
                        .percent(signature.percent())
                        .build())
                .collect(Collectors.toList());
    }

    @NotNull
    public static String annotateWithEtiology(@NotNull com.hartwig.hmftools.common.sigs.SignatureAllocation signatureAllocation,
            @NotNull Map<String, String> etiologyPerSignature)
    {
        String signatureName = signatureAllocation.signature();
        if(etiologyPerSignature.containsKey(signatureName))
        {
            return etiologyPerSignature.get(signatureName);
        }
        else if("MISALLOC".equals(signatureName))
        {
            return "-";
        }
        else
        {
            throw new IllegalArgumentException("Signature " + signatureName + " is not in signatures etiology tsv");
        }
    }
}
