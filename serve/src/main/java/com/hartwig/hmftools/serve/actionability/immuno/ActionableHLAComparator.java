package com.hartwig.hmftools.serve.actionability.immuno;

import java.util.Comparator;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.ActionableEventComparator;

import org.jetbrains.annotations.NotNull;

class ActionableHLAComparator {

}
//class ActionableHLAComparator implements Comparator<ActionableHLAComparator> {
//
//    @NotNull
//    private final Comparator<ActionableEvent> actionableEventComparator = new ActionableEventComparator();
//
//    @Override
//    public int compare(@NotNull ActionableHLA hla1, @NotNull ActionableHLA hla2) {
//        int hlaCompare = hla1.hlaTypering().compareTo(hla2.hlaTypering());
//        if (hlaCompare != 0) {
//            return hlaCompare;
//        }
//
//        return actionableEventComparator.compare(hla1, hla2);
//    }
//}
