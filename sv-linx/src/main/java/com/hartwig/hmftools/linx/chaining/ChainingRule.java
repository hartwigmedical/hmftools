package com.hartwig.hmftools.linx.chaining;

import java.util.List;
import static java.lang.Math.pow;

public enum ChainingRule
{
    ASSEMBLY,
    SINGLE_OPTION,
    CHAIN_SPLIT,
    FOLDBACK,
    PLOIDY_MATCH,
    PLOIDY_OVERLAP,
    ADJACENT,
    PLOIDY_MAX,
    NEAREST;

    public static int ruleToPriority(ChainingRule rule)
    {
        switch(rule)
        {
            case SINGLE_OPTION: return 8;
            case CHAIN_SPLIT: return 7;
            case FOLDBACK: return 6;
            case PLOIDY_MATCH: return 5;
            case PLOIDY_OVERLAP: return 4;
            case ADJACENT: return 3;
            case PLOIDY_MAX: return 2;
            case NEAREST: return 1;
            default: return 0;
        }
    }

    public static int calcRulePriority(final List<ChainingRule> rules)
    {
        return rules.stream().mapToInt(x -> (int)pow(2, ruleToPriority(x))).sum();
    }
}
