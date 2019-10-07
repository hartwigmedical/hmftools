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
    public abstract String targeting();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract List<String> cDNA();

    @NotNull
    public abstract String primary_tumor_type();

    @NotNull
    public abstract List<String> individual_mutation();

    @NotNull
    public abstract String drugsFullName();

    @NotNull
    public abstract String curator();

    @NotNull
    public abstract String drug_family();

    @NotNull
    public abstract String alteration();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String biomarker();

    @NotNull
    public abstract List<String> gDNA();

    @NotNull
    public abstract String drug_status();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract List<String> transcript();

    @NotNull
    public abstract List<String> strand();

    @NotNull
    public abstract List<String> info();

    @NotNull
    public abstract String assay_type();

    @NotNull
    public abstract String alteration_type();

    @NotNull
    public abstract List<String> region();

    @NotNull
    public abstract String evidence_level();

    @NotNull
    public abstract String association();

    @NotNull
    public abstract String metastatic_Tumor_Type();

}
