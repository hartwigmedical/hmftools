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
        public abstract String type();

        public abstract String start();

        public abstract String geneStart();

        public abstract String geneContextStart();

        public abstract String transcriptStart();

        public abstract String end();

        public abstract String geneEnd();

        public abstract String geneContextEnd();

        public abstract String transcriptEnd();

        public abstract String vaf();
    }

    @Value.Immutable
    public static abstract class GeneDisruption {
        public abstract String geneName();

        public abstract String location();

        public abstract String geneContext();

        public abstract String transcript();

        public abstract String partner();

        public abstract String type();

        public abstract String orientation();

        public abstract String vaf();
    }

    @NotNull
    public abstract List<VariantAnnotation> annotations();

    @NotNull
    public abstract List<GeneFusion> fusions();

    @NotNull
    public abstract List<GeneDisruption> disruptions();
}
