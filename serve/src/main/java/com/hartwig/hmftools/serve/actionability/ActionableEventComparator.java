package com.hartwig.hmftools.serve.actionability;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

public class ActionableEventComparator implements Comparator<ActionableEvent> {

    @Override
    public int compare(@NotNull ActionableEvent event1, @NotNull ActionableEvent event2) {
        int sourceCompare = event1.source().toString().compareTo(event2.source().toString());
        if (sourceCompare != 0) {
            return sourceCompare;
        }

        int levelCompare = event1.level().toString().compareTo(event2.level().toString());
        if (levelCompare != 0) {
            return levelCompare;
        }

        int directionCompare = event1.direction().toString().compareTo(event2.direction().toString());
        if (directionCompare != 0) {
            return directionCompare;
        }

        int treatmentCompare = event1.treatment().compareTo(event2.treatment());
        if (treatmentCompare != 0) {
            return treatmentCompare;
        }

        int cancerTypeCompare = event1.applicableCancerType().cancerType().compareTo(event2.applicableCancerType().cancerType());
        if (cancerTypeCompare != 0) {
            return cancerTypeCompare;
        }

        return event1.applicableCancerType().doid().compareTo(event2.applicableCancerType().doid());
    }
}