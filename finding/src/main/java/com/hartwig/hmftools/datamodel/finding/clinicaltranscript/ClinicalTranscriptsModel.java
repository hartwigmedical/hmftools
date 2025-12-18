package com.hartwig.hmftools.datamodel.finding.clinicaltranscript;

import java.util.Map;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ClinicalTranscriptsModel {

    @NotNull
    private final Map<String, String> clinicalTranscriptMap;

    public static ClinicalTranscriptsModel emptyModel() {
        return new ClinicalTranscriptsModel(Map.of());
    }

    ClinicalTranscriptsModel(@NotNull final Map<String, String> clinicalTranscriptMap) {
        this.clinicalTranscriptMap = clinicalTranscriptMap;
    }

    @Nullable
    public String findCanonicalTranscriptForGene(@NotNull String gene) {
        return clinicalTranscriptMap.get(gene);
    }
}
