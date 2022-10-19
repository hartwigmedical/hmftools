package com.hartwig.hmftools.common.serve.actionability.immuno;

import java.util.Comparator;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.actionability.ActionableEventComparator;

import org.jetbrains.annotations.NotNull;

class ActionableHLAComparator implements Comparator<ActionableHLA> {

    @NotNull
    private final Comparator<ActionableEvent> actionableEventComparator = new ActionableEventComparator();

    @Override
    public int compare(@NotNull ActionableHLA hla1, @NotNull ActionableHLA hla2) {
        int hlaCompare = hla1.hlaType().compareTo(hla2.hlaType());
        if (hlaCompare != 0) {
            return hlaCompare;
        }

        return actionableEventComparator.compare(hla1, hla2);
    }
}