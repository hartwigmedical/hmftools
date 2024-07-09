package com.hartwig.hmftools.datamodel.linx;

import org.jetbrains.annotations.NotNull;

public enum LinxReportableReason
{
    OK(""),
    NOT_KNOWN("Unknown fusion pair"),
    UNPHASED_NOT_KNOWN("Unphased, no known fusion pair"),
    UNPHASED_5P_UTR("Unphased, 5P UTR"),
    UNPHASED_SHORT("Unphased, short unphased distance"),
    SGL_NOT_KNOWN("SGL, no known fusion pair"),
    PRE_GENE_DISTANCE("Max upstream distance"),
    NMD("Nonsense mediated decay"),
    NEG_SPLICE_ACC_DISTANCE("Negative previous splice acceptor distance"),
    EXON_SKIPPING("Exon skipping"),
    CHAIN_TERMINATED("Chain terminated"),
    NON_DISRUPTIVE_CHAIN("Non-disruptive chain"),
    INVALID_TRAVERSAL("Invalid chain traversal"),
    CHAIN_LINKS("Maximum chain links"),
    PROTEIN_DOMAINS("Protein domains");

    private final @NotNull String display;

    LinxReportableReason(@NotNull String display)
    {
        this.display = display;
    }

    @NotNull
    public String display()
    {
        return display;
    }
}
