package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class DisruptionAnalyzer {

    @NotNull
    private final List<DisruptionReaderFile> disruptionReaderFiles;

    DisruptionAnalyzer(@NotNull final List<DisruptionReaderFile> disruptionReaderFiles) {
        this.disruptionReaderFiles = disruptionReaderFiles;
    }

    @NotNull
    public List<DisruptionReaderFile> filteringDisruptions() {
        List<DisruptionReaderFile> raportableDisruptions = Lists.newArrayList();
        for (DisruptionReaderFile disruption: disruptionReaderFiles) {
            if (disruption.reportable().equals(true)) {
                raportableDisruptions.add(disruption);
            }
        }
        return raportableDisruptions;
    }
}
