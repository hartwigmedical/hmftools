package com.hartwig.hmftools.common.serve.actionability.characteristic;

import java.util.Comparator;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.actionability.ActionableEventComparator;

import org.jetbrains.annotations.NotNull;

class ActionableCharacteristicComparator implements Comparator<ActionableCharacteristic> {

    @NotNull
    private final Comparator<ActionableEvent> actionableEventComparator = new ActionableEventComparator();

    @Override
    public int compare(@NotNull ActionableCharacteristic characteristic1, @NotNull ActionableCharacteristic characteristic2) {
        int nameCompare = characteristic1.name().toString().compareTo(characteristic2.name().toString());
        if (nameCompare != 0) {
            return nameCompare;
        }

        return actionableEventComparator.compare(characteristic1, characteristic2);
    }
}
