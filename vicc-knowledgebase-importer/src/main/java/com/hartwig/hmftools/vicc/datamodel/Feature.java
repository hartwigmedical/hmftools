package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import com.hartwig.hmftools.vicc.util.EventAnnotationExtractor;
import com.hartwig.hmftools.vicc.util.ProteinAnnotationExtractor;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Feature {

    @NotNull
    @Value.Derived
    public String proteinAnnotation() {
        return ProteinAnnotationExtractor.toProteinAnnotation(name());
    }

    @NotNull
    @Value.Derived
    public String eventAnnotation() {
        return EventAnnotationExtractor.toEventAnnotation();
    }

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String biomarkerType();

    @Nullable
    public abstract String referenceName();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String start();

    @Nullable
    public abstract String end();

    @Nullable
    public abstract String ref();

    @Nullable
    public abstract String alt();

    @NotNull
    public abstract List<String> provenance();

    @Nullable
    public abstract String provenanceRule();

    @Nullable
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String entrezId();

    @Nullable
    public abstract SequenceOntology sequenceOntology();

    @NotNull
    public abstract List<String> links();

    @Nullable
    public abstract String description();

    @Nullable
    public abstract FeatureInfo info();

    @Nullable
    public abstract FeatureAttribute attribute();
}
