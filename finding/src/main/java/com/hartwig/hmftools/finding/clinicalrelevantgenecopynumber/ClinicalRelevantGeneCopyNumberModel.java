package com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class ClinicalRelevantGeneCopyNumberModel
{
    @NotNull
    private final List<String> clinicalRelevantGeneCopyNumbersList;

    ClinicalRelevantGeneCopyNumberModel(@NotNull final List<String> clinicalRelevantGeneCopyNumbersList) {
        this.clinicalRelevantGeneCopyNumbersList = clinicalRelevantGeneCopyNumbersList;
    }

    public boolean findClinicalRelevantCopyNumberGene(@NotNull String gene) {
        int indexGene = clinicalRelevantGeneCopyNumbersList.indexOf(gene);
        return indexGene != -1;
    }
}