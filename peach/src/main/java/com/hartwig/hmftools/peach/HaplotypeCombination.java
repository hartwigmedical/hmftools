package com.hartwig.hmftools.peach;

import org.jetbrains.annotations.NotNull;

import java.util.Map;

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;

public class HaplotypeCombination
{
    @NotNull
    Map<String, Integer> haplotypeNameToCount;

    public HaplotypeCombination(@NotNull Map<String, Integer> haplotypeNameToCount)
    {
        if (haplotypeNameToCount.values().stream().anyMatch(c -> c <= 0))
        {
            String error_msg = String.format(
                    "Invalid haplotype combination with non-positive count(s): %s", haplotypeNameToCount
            );
            throw new RuntimeException(error_msg);
        }
        else if (haplotypeNameToCount.values().stream().mapToInt(i -> i).sum() < GERMLINE_TOTAL_COPY_NUMBER)
        {
            String error_msg = String.format(
                    "Less than %d haplotype calls in combination is not allowed: %s",
                    GERMLINE_TOTAL_COPY_NUMBER, haplotypeNameToCount
            );
            throw new RuntimeException(error_msg);
        }
        this.haplotypeNameToCount = haplotypeNameToCount;
    }

    public String toString()
    {
        return String.format("HaplotypeCombination(%s)", haplotypeNameToCount);
    }
}
