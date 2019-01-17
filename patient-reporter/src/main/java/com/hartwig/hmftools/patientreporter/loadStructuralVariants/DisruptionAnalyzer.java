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
        for (int i = 1; i < disruptionReaderFiles.size(); i ++) {
            if (disruptionReaderFiles.get(i).reportable().equals(true)) {
                raportableDisruptions.add(disruptionReaderFiles.get(i));
            }
        }
        return raportableDisruptions;
    }
}
