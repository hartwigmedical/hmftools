package com.hartwig.hmftools.common.serve.actionability.range;

import java.util.Comparator;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.actionability.ActionableEventComparator;

import org.jetbrains.annotations.NotNull;

class ActionableRangeComparator implements Comparator<ActionableRange> {

    @NotNull
    private final Comparator<ActionableEvent> actionableEventComparator = new ActionableEventComparator();

    @Override
    public int compare(@NotNull ActionableRange range1, @NotNull ActionableRange range2) {
        int defaultCompare = range1.compareTo(range2);
        if (defaultCompare != 0) {
            return defaultCompare;
        }

        int geneCompare = range1.gene().compareTo(range2.gene());
        if (geneCompare != 0) {
            return geneCompare;
        }

        return actionableEventComparator.compare(range1, range2);
    }
}
