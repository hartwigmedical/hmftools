package com.hartwig.hmftools.bachelor;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class EligibilityReport {
    @NotNull
    public abstract String patient();

    enum ReportType {
        GERMLINE_MUTATION,
        SOMATIC_MUTATION,
        GERMLINE_DELETION,
        SOMATIC_DELETION,
        SOMATIC_DISRUPTION
    }

    enum MatchType {
        NONE,
        GENE_TRANSCRIPT,
        HOTSPOT_LOCATION,
        WHITELIST
    }

    @NotNull
    public abstract ReportType source();

    public abstract MatchType matchType();

    @NotNull
    public abstract String program();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String genes();

    @NotNull
    public abstract String transcriptId();

    @NotNull
    public abstract String chrom();

    public abstract long pos();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alts();

    @NotNull
    public abstract String effects();

    @NotNull
    public abstract String annotations();

    public abstract String hgvsProtein();

    public abstract String hgvsCoding();

    public abstract boolean isHomozygous();

    public abstract int phredScore();

    public abstract int altCount();

    public abstract int readDepth();
}