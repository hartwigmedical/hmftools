package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DisruptionAnalyzer {

    @NotNull
    private final List<DisruptionReaderFile> disruptionReaderFiles;

    DisruptionAnalyzer(@NotNull final List<DisruptionReaderFile> disruptionReaderFiles) {
        this.disruptionReaderFiles = disruptionReaderFiles;
    }

    @Nullable
    public List<DisruptionReaderFile> filteringDisruptions() {
        List<DisruptionReaderFile> raportableDisruptions = Lists.newArrayList();
        for (int i = 1; i <= disruptionReaderFiles.size(); i ++) {
        }
        return raportableDisruptions;
    }
}
