package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.CANDIDATE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.NO_HET_LOCI;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.SOLUTION;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNMATCHED_AMINO_ACID;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.WILD_ONLY;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CoverageQC {

    public final int ATypes;
    public final int BTypes;
    public final int CTypes;

    public final int TotalFragments; // fragment with at least 1 sufficient map-qual base

    public final int FittedFragments; // fragments allocated to winning solution
    public final int UnmatchedFragments;
    public final int UninformativeFragments;
    public final int HlaYFragments;


    // winning solution breakdown
    public final int UniqueFragments; //
    public final int SharedFragments;
    public final int WildcardFragments;

    public final double PercentUnique;
    public final double PercentShared;
    public final double PercentWildcard;

    /*
    - total fragments - of fragments overlapping a coding base of HLA-A, HLA-B or HLA-C with sufficient MAPQ
    - FittedFragments - total fragments allocated exactly to solution (ie. exact or wildcard match at every heterozygous location)
    - discardedHLA_YFragments - fragments matching HLA-Y
    - uninformativeFragments - fragment does not overlap any heterozygous location considered in evidence phase
    - unmatchedFragments - fragment does not match any allele exactly in solution at all heterozygous location (eg. due to sequencing error or mapping error)
    - aTypes/bTypes/cTypes - # of distinct HLA-A alleles fitted (0,1 or 2), etc
    - proportionUnique - Percentage of fitted fragments that are uniquely assigned to 1 allele
    - proportionShared - Percentage of fitted fragments allocated across multiple alleles
    - proportionWildcard - Percentage of fitted fragments uniquely assigned to wildcard regions
    - discardedAlignmentFragments - fragments discarded because 1 read aligns more than 1000 bases from a HLA A,B or C gene
    - discardedIndelFragments - count of fragments excluded due to not matching a known indel
    - discardedIndelMaxSupport - Maximum fragment support for an indel detected but not present in any known allele
     */


    public CoverageQC(
            int aTypes, int bTypes, int cTypes, int totalFragments, int unmatched, int uninformative, int hlaY,
            int uniqueFragments, int sharedFragments, int wildcardFragments)
    {
        ATypes = aTypes;
        BTypes = bTypes;
        CTypes = cTypes;
        TotalFragments = totalFragments;
        UnmatchedFragments = unmatched;
        UninformativeFragments = uninformative;
        HlaYFragments = hlaY;

        UniqueFragments = uniqueFragments;
        SharedFragments = sharedFragments;
        WildcardFragments = wildcardFragments;

        FittedFragments = uniqueFragments + sharedFragments + wildcardFragments;
        PercentUnique = 1.0 * uniqueFragments / FittedFragments;
        PercentShared = 1.0 * sharedFragments / FittedFragments;
        PercentWildcard = 1.0 * wildcardFragments / FittedFragments;
    }

    public static CoverageQC create(
            final List<Fragment> fragments, final HlaComplexCoverage winner)
    {
        List<HlaAllele> alleles = winner.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList());
        int aTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_A)).collect(Collectors.toSet()).size();
        int bTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_B)).collect(Collectors.toSet()).size();
        int cTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_C)).collect(Collectors.toSet()).size();

        if (aTypes == 0 || bTypes == 0 || cTypes == 0)
        {
            LL_LOGGER.warn("  UNMATCHED ALLELE: {} A alleles, {} B alleles, {} C alleles", aTypes, bTypes, cTypes);
        }

        /*
        - total fragments - of fragments overlapping a coding base of HLA-A, HLA-B or HLA-C with sufficient MAPQ
        - FittedFragments - total fragments allocated exactly to solution (ie. exact or wildcard match at every heterozygous location)
        - discardedHLA_YFragments - fragments matching HLA-Y
        - uninformativeFragments - fragment does not overlap any heterozygous location considered in evidence phase
        - unmatchedFragments - fragment does not match any allele exactly in solution at all heterozygous location (eg. due to sequencing error or mapping error)
        - aTypes/bTypes/cTypes - # of distinct HLA-A alleles fitted (0,1 or 2), etc
        - proportionUnique - Percentage of fitted fragments that are uniquely assigned to 1 allele
        - proportionShared - Percentage of fitted fragments allocated across multiple alleles
        - proportionWildcard - Percentage of fitted fragments uniquely assigned to wildcard regions
        - discardedAlignmentFragments - fragments discarded because 1 read aligns more than 1000 bases from a HLA A,B or C gene
        - discardedIndelFragments - count of fragments excluded due to not matching a known indel
        - discardedIndelMaxSupport - Maximum fragment support for an indel detected but not present in any known allele
         */

        int totalFragments = fragments.size();

        int solutionFragments = 0; // winner.TotalCoverage;
        int hlaYFragments = 0;
        int uninformativeFragments = 0;
        int unmatchedFragments = 0;

        for(Fragment fragment : fragments)
        {
            if(fragment.scope() == HLA_Y)
                ++hlaYFragments;
            else if(fragment.scope() == SOLUTION)
                ++solutionFragments;
            else if(fragment.scope() == NO_HET_LOCI)
                ++uninformativeFragments;
            else if(fragment.scope().isUnmatched())
                ++unmatchedFragments;
            else
                LL_LOGGER.warn("frag({} : {}) scope({}) unassigned", fragment.id(), fragment.readInfo(), fragment.scope());
        }

        int unassignedFragments = totalFragments - solutionFragments - hlaYFragments - uninformativeFragments - unmatchedFragments;

        if(unassignedFragments != 0)
        {
            LL_LOGGER.warn("fragment count inconsistency: total({}) sol({}) unmatch({}) uninform({}) hlaY({}) unassigned({})",
                    totalFragments, solutionFragments, unmatchedFragments, uninformativeFragments, hlaYFragments, unassignedFragments);
        }

        return new CoverageQC(
                aTypes, bTypes, cTypes,
                totalFragments, unmatchedFragments, uninformativeFragments, hlaYFragments,
                winner.UniqueCoverage, winner.SharedCoverage, winner.WildCoverage);
    }

    public List<String> header()
    {
        return Lists.newArrayList(
                "ATypes", "BTypes", "CTypes",
                "TotalFragments", "FittedFragments", "UnmatchedFragments", "UninformativeFragments", "HlaYFragments",
                "PercentUnique", "PercentShared", "PercentWildcard");
    }

    public List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(ATypes), String.valueOf(BTypes), String.valueOf(CTypes),
                String.valueOf(TotalFragments), String.valueOf(FittedFragments),
                String.valueOf(UnmatchedFragments), String.valueOf(UninformativeFragments), String.valueOf(HlaYFragments),
                String.format("%.3f", PercentUnique), String.format("%.3f", PercentShared), String.format("%.3f", PercentWildcard));
    }

    /*
    public String toString()
    {
        return String.format("CoverageQC(aTypes=%d, bTypes=%d, cTypes=%d, unusedFragments=%d, uniqueFragments=%d, sharedFragments=%d, "
                        + "wildcardFragments=%d, fittedFragments=%d, percentUnique=%.3f, percentShared=%.3f, percentWildcard=%.3f)",
                ATypes, BTypes, CTypes, UnmatchedFragments, UniqueFragments, SharedFragments, FittedFragments,
                PercentUnique, PercentShared, PercentWildcard);
    }
     */
}