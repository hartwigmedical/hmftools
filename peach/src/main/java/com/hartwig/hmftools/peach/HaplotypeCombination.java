package com.hartwig.hmftools.peach;

import org.jetbrains.annotations.NotNull;

import java.util.Map;

public class HaplotypeCombination
{
    @NotNull
    Map<String, Integer> haplotypeNameToCount;

    public HaplotypeCombination(@NotNull Map<String, Integer> haplotypeNameToCount)
    {
        if (haplotypeNameToCount.values().stream().anyMatch(c -> c < 0))
        {
            String error_msg = String.format(
                    "Invalid haplotype combination with negative count(s): %s", haplotypeNameToCount
            );
            throw new RuntimeException(error_msg);
        }
        else if (haplotypeNameToCount.values().stream().mapToInt(i -> i).sum() < 2)
        {
            String error_msg = String.format(
                    "Less than 2 haplotype calls in combination is not allowed: %s", haplotypeNameToCount
            );
            throw new RuntimeException(error_msg);
        }
        this.haplotypeNameToCount = haplotypeNameToCount;
    }
}
