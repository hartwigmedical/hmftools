package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record Fusion(
        @NotNull DriverFields driver,
        @NotNull String geneStart,
        @NotNull String geneContextStart,
        @NotNull String geneTranscriptStart,
        @NotNull String geneEnd,
        @NotNull String geneContextEnd,
        @NotNull String geneTranscriptEnd,
        @NotNull FusionType reportedType,
        @NotNull List<UnreportableReason> unreportedReasons,
        @NotNull FusionPhasedType phased,
        int fusedExonUp,
        int fusedExonDown,
        int chainLinks,
        boolean chainTerminated,
        @NotNull List<String> domainsKept,
        @NotNull List<String> domainsLost,
        double junctionCopyNumber
) implements Driver
{
    public enum FusionType
    {
        NONE,
        PROMISCUOUS_3,
        PROMISCUOUS_5,
        PROMISCUOUS_BOTH,
        IG_PROMISCUOUS,
        KNOWN_PAIR,
        IG_KNOWN_PAIR,
        EXON_DEL_DUP,
        PROMISCUOUS_ENHANCER_TARGET
    }

    public enum UnreportableReason
    {
        NONE,
        NOT_KNOWN,
        UNPHASED_NOT_KNOWN,
        UNPHASED_5P_UTR,
        UNPHASED_SHORT,
        SGL_NOT_KNOWN,
        PRE_GENE_DISTANCE,
        NONSENSE_MEDIATED_DECAY,
        NEG_SPLICE_ACC_DISTANCE,
        EXON_SKIPPING,
        CHAIN_TERMINATED,
        NON_DISRUPTIVE_CHAIN,
        INVALID_TRAVERSAL,
        CHAIN_LINKS,
        DISRUPTED_PROTEIN_DOMAINS
    }

    public enum FusionPhasedType
    {
        INFRAME,
        SKIPPED_EXONS,
        OUT_OF_FRAME
    }

    @NotNull
    @Override
    public String findingKey()
    {
        return driver.findingKey();
    }

    @NotNull
    @Override
    public DriverSource driverSource()
    {
        return driver.driverSource();
    }

    @NotNull
    @Override
    public ReportedStatus reportedStatus()
    {
        return driver.reportedStatus();
    }

    @NotNull
    @Override
    public DriverInterpretation driverInterpretation()
    {
        return driver.driverInterpretation();
    }

    @Override
    public double driverLikelihood()
    {
        return driver.driverLikelihood();
    }

    @NotNull
    public String display()
    {
        return String.format("%s::%s", geneStart, geneEnd);
    }
}
