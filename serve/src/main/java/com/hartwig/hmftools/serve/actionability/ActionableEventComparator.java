package com.hartwig.hmftools.serve.actionability;

import java.util.Comparator;

import com.hartwig.hmftools.serve.cancertype.CancerType;

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

        return compareCancerTypes(event1.applicableCancerType(), event2.applicableCancerType());
    }

    private static int compareCancerTypes(@NotNull CancerType cancerType1, @NotNull CancerType cancerType2) {
        int nameCompare = cancerType1.name().compareTo(cancerType2.name());
        if (nameCompare != 0) {
            return nameCompare;
        }

        return cancerType1.doid().compareTo(cancerType2.doid());
    }
}