package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Feature {

    @NotNull
    @Value.Derived
    public String proteinAnnotation() {
        String trimmedName = name().trim();
        // Many KBs include the gene in the feature name in some form (eg "EGFR E709K" or "EGFR:E709K").
        // Other KBs put the coding info behind the protein annotation ("V130L (c.388G>C)" rather than the gene in front of it)
        String proteinAnnotation;
        if (trimmedName.contains(" ")) {
            String[] trimmedParts = trimmedName.split(" ");
            proteinAnnotation = trimmedParts[1].contains("(c.") ? trimmedParts[0] : trimmedParts[1];
        } else if (trimmedName.contains(":")) {
            proteinAnnotation = trimmedName.split(":")[1];
        } else {
            proteinAnnotation = trimmedName;
        }

        // Some KBs include "p." in front of the protein annotation
        return proteinAnnotation.startsWith("p.") ? proteinAnnotation.substring(2) : proteinAnnotation;
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
