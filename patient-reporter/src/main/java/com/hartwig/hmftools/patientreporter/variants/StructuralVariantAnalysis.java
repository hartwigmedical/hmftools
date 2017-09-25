package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.svannotation.VariantAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantAnalysis {

    @Value.Immutable
    public static abstract class GeneFusion {
        public abstract String Type();
        public abstract String Start();
        public abstract String GeneStart();
        public abstract String GeneContextStart();
        public abstract String TranscriptStart();
        public abstract String End();
        public abstract String GeneEnd();
        public abstract String GeneContextEnd();
        public abstract String TranscriptEnd();
        public abstract String VAF();
    }

    @Value.Immutable
    public static abstract class GeneDisruption {
        public abstract String GeneName();
        public abstract String Location();
        public abstract String GeneContext();
        public abstract String Transcript();
        public abstract String Partner();
        public abstract String Type();
        public abstract String Orientation();
        public abstract String VAF();
    }

    @NotNull
    public abstract List<VariantAnnotation> annotations();

    @NotNull
    public abstract List<GeneFusion> fusions();

    @NotNull
    public abstract List<GeneDisruption> disruptions();
}
