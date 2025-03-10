package com.hartwig.hmftools.peach.haplotype;

import java.util.HashMap;
import java.util.Map;

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;

public class HaplotypeCombination
{
    private final Map<String, Integer> haplotypeNameToCount;

    public HaplotypeCombination(final Map<String, Integer> haplotypeNameToCount)
    {
        if(haplotypeNameToCount.values().stream().anyMatch(c -> c <= 0))
        {
            String errorMsg = String.format("Invalid haplotype combination with non-positive count(s): %s", haplotypeNameToCount);
            throw new RuntimeException(errorMsg);
        }
        else if(haplotypeNameToCount.values().stream().mapToInt(i -> i).sum() < GERMLINE_TOTAL_COPY_NUMBER)
        {
            String errorMsg = String.format(
                    "Less than %d haplotype calls in combination is not allowed: %s",
                    GERMLINE_TOTAL_COPY_NUMBER, haplotypeNameToCount
            );
            throw new RuntimeException(errorMsg);
        }
        this.haplotypeNameToCount = haplotypeNameToCount;
    }

    public boolean equals(final Object other)
    {
        if(this == other)
        {
            return true;
        }

        if(!(other instanceof HaplotypeCombination))
        {
            return false;
        }

        return haplotypeNameToCount.equals(((HaplotypeCombination) other).haplotypeNameToCount);
    }

    public String toString()
    {
        return String.format("HaplotypeCombination(%s)", haplotypeNameToCount);
    }

    public Map<String, Integer> getHaplotypeNameToCount()
    {
        return new HashMap<>(haplotypeNameToCount);
    }

    public int getHaplotypeCount()
    {
        return haplotypeNameToCount.values().stream().mapToInt(c -> c).sum();
    }

    public int getHaplotypeCountWithout(final String ignoredHaplotypeName)
    {
        return haplotypeNameToCount.entrySet()
                .stream()
                .filter(e -> !e.getKey().equals(ignoredHaplotypeName))
                .mapToInt(Map.Entry::getValue)
                .sum();
    }
}
