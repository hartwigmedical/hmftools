package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvAnalyzerModel {


    @NotNull
    private final DisruptionAnalyzer disruptionAnalyzer;
    @NotNull
    private final FusionAnalyzer fusionAnalyzer;

    @NotNull
    public static SvAnalyzerModel fromKnowledgebase(@NotNull String fusionFile, @NotNull String disruptionFile) throws IOException {
        DisruptionAnalyzer disruptionAnalyzer = DisruptionFactory.readingDisruption(disruptionFile);
        FusionAnalyzer fusionAnalyzer = FusionFactory.readingFusion(fusionFile);
        return new SvAnalyzerModel(fusionAnalyzer, disruptionAnalyzer);
    }

    private SvAnalyzerModel(@NotNull final FusionAnalyzer fusionAnalyzer,
            @NotNull final DisruptionAnalyzer disruptionAnalyzer) {
        this.fusionAnalyzer = fusionAnalyzer;
        this.disruptionAnalyzer = disruptionAnalyzer;
    }

    @Nullable
    public List<FusionReaderFile> filterFusions() {
        return fusionAnalyzer.filteringFusions();

    }

    @Nullable
    public List<DisruptionReaderFile> filterDisruptions() {
        return disruptionAnalyzer.filteringDisruptions();
    }
}
