package com.hartwig.hmftools.serve.actionability.characteristic;

import java.util.Comparator;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.ActionableEventComparator;

import org.jetbrains.annotations.NotNull;

class ActionableSignatureComparator implements Comparator<ActionableSignature> {

    @NotNull
    private final Comparator<ActionableEvent> actionableEventComparator = new ActionableEventComparator();

    @Override
    public int compare(@NotNull ActionableSignature signature1, @NotNull ActionableSignature signature2) {
        int signatureCompare = signature1.name().toString().compareTo(signature2.name().toString());
        if (signatureCompare != 0) {
            return signatureCompare;
        }

        return actionableEventComparator.compare(signature1, signature2);
    }
}
