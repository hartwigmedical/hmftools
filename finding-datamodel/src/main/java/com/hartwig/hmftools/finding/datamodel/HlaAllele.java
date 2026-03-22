package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

// HLA allele nomenclature:
//
//   HLA-A  *  02  : 101
//      |_|    |_|   |_|
//       |      |     |
//   geneSymbol |   hlaProtein
//          alleleGroup
//
// geneSymbol  - e.g. "A", "B", "DQA1"
// alleleGroup - two-digit field identifying the allele group (serological antigen), e.g. "02"
// hlaProtein  - two or three digit field identifying the specific protein, e.g. "101"
@RecordBuilder
public record HlaAllele(
        @NotNull String findingKey,
        @NotNull GeneClass geneClass,
        @NotNull String geneSymbol,
        @NotNull String alleleGroup,
        @NotNull String hlaProtein,
        @NotNull Set<QcStatus> qcStatus,
        int germlineCopyNumber,
        @Nullable Double tumorCopyNumber,
        @Nullable Integer refFragments,
        int tumorFragments,
        @Nullable Integer rnaFragments,
        double somaticMissense,
        double somaticNonsenseOrFrameshift,
        double somaticSplice,
        double somaticSynonymous,
        double somaticInframeIndel
) implements Finding
{
    public enum GeneClass
    {
        MHC_CLASS_I,
        MHC_CLASS_II
    }

    public enum QcStatus
    {
        PASS,
        FAIL,
        WARN_UNMATCHED_ALLELE,
        WARN_UNMATCHED_SOMATIC_VARIANT,
        WARN_UNMATCHED_HAPLOTYPE,
        WARN_UNMATCHED_AMINO_ACID,
        WARN_LOW_COVERAGE,
        WARN_LOW_BASE_QUAL,
        WARN_UNMATCHED_INDEL
    }

    public @NotNull String gene()
    {
        return "HLA-" + geneSymbol;
    }

    public @NotNull String allele()
    {
        return geneSymbol + '*' + alleleGroup + ':' + hlaProtein;
    }

    public boolean hasSomaticVariants()
    {
        return Doubles.positive(somaticMissense()) || Doubles.positive(somaticNonsenseOrFrameshift()) || Doubles.positive(
                somaticSplice()) || Doubles.positive(somaticInframeIndel());
    }
}
