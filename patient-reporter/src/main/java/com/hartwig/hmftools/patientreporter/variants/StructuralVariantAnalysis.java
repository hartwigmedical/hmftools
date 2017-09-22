package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.svannotation.VariantAnnotation;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalysis {

    public static class GeneFusion {
        public String Type;
        public String Start;
        public String GeneStart;
        public String GeneContextStart;
        public String TranscriptStart;
        public String End;
        public String GeneEnd;
        public String GeneContextEnd;
        public String TranscriptEnd;
        public String VAF;
    }

    public static class GeneDisruption {
        public String GeneName;
        public String Location;
        public String GeneContext;
        public String Transcript;
        public String Partner;
        public String HGVS;
        public String Type;
        public String Orientation;
        public String VAF;
        public String TAF = "TODO";
    }

    @NotNull
    private final List<VariantAnnotation> annotations;
    @NotNull
    private final List<GeneFusion> fusions;
    @NotNull
    private final List<GeneDisruption> disruptions;

    StructuralVariantAnalysis(@NotNull final List<VariantAnnotation> annotations, @NotNull final List<GeneFusion> fusions,
            @NotNull final List<GeneDisruption> disruptions) {
        this.annotations = annotations;
        this.fusions = fusions;
        this.disruptions = disruptions;
    }

    @NotNull
    public List<VariantAnnotation> annotations() {
        return ImmutableList.copyOf(annotations);
    }

    @NotNull
    public List<GeneFusion> fusions() {
        return ImmutableList.copyOf(fusions);
    }

    @NotNull
    public List<GeneDisruption> disruptions() {
        return ImmutableList.copyOf(disruptions);
    }
}
