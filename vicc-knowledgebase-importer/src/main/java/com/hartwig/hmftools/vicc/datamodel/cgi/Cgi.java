package com.hartwig.hmftools.vicc.datamodel.cgi;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Cgi implements KbSpecificObject {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String biomarker();

    @NotNull
    public abstract String alteration();

    @NotNull
    public abstract String alterationType();

    @NotNull
    public abstract List<String> transcripts();

    @NotNull
    public abstract List<String> individualMutations();

    @NotNull
    public abstract List<String> gDNA();

    @NotNull
    public abstract List<String> cDNA();

    @NotNull
    public abstract List<String> info();

    @NotNull
    public abstract List<String> regions();

    @NotNull
    public abstract List<String> strands();

    @NotNull
    public abstract String association();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String drugFamily();

    @NotNull
    public abstract String drugFullName();

    @NotNull
    public abstract String drugStatus();

    @NotNull
    public abstract String targeting();

    @NotNull
    public abstract String primaryTumorType();

    @NotNull
    public abstract String metastaticTumorType();

    @NotNull
    public abstract String evidenceLevel();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String curator();

    @NotNull
    public abstract String assayType();

}
