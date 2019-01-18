package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvAnalyzerModel {

    @NotNull
    private final FusionAnalyzer fusionAnalyzer;
    @NotNull
    private final DisruptionAnalyzer disruptionAnalyzer;

    @NotNull
    public static SvAnalyzerModel fromFiles(@NotNull String fusionFile, @NotNull String disruptionFile) throws IOException {
        FusionAnalyzer fusionAnalyzer = FusionFactory.fromFusionFile(fusionFile);
        DisruptionAnalyzer disruptionAnalyzer = DisruptionFactory.fromDisruptionFile(disruptionFile);
        return new SvAnalyzerModel(fusionAnalyzer, disruptionAnalyzer);
    }

    private SvAnalyzerModel(@NotNull final FusionAnalyzer fusionAnalyzer, @NotNull final DisruptionAnalyzer disruptionAnalyzer) {
        this.fusionAnalyzer = fusionAnalyzer;
        this.disruptionAnalyzer = disruptionAnalyzer;
    }

    @Nullable
    public List<Fusion> filterFusions() {
        return fusionAnalyzer.reportableFusions();
    }

    @Nullable
    public List<Disruption> filterDisruptions(@NotNull GeneModel geneModel) {
        return disruptionAnalyzer.reportableDisruptions(geneModel);
    }
}
