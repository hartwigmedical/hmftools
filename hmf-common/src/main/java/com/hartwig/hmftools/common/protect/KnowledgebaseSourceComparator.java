package com.hartwig.hmftools.common.protect;

import org.jetbrains.annotations.NotNull;

import java.util.Comparator;

public class KnowledgebaseSourceComparator implements Comparator<KnowledgebaseSource> {

    @Override
    public int compare(@NotNull KnowledgebaseSource source1, @NotNull KnowledgebaseSource source2) {

        int knowledgebaseNameCompare = source1.name().toString().compareTo(source2.name().toString());
        if (knowledgebaseNameCompare != 0) {
            return knowledgebaseNameCompare;
        }

        int sourceEventCompare = source1.sourceEvent().compareTo(source2.sourceEvent());
        if (sourceEventCompare != 0) {
            return sourceEventCompare;
        }

        int evidenceTypeCompare = source1.evidenceType().display().compareTo(source2.evidenceType().display());
        if (evidenceTypeCompare != 0) {
            return evidenceTypeCompare;
        }

        return 0;
    }
}