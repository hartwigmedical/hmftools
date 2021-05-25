package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CoverageQC {

    public final int ATypes;
    public final int BTypes;
    public final int CTypes;
    public final int TotalFragments;
    public final int UniqueFragments;
    public final int SharedFragments;
    public final int WildcardFragments;

    public final int FittedFragments;
    public final int UnusedFragments;
    public final double PercentUnique;
    public final double PercentShared;
    public final double PercentWildcard;

    public CoverageQC(
            final int aTypes, final int bTypes, final int cTypes, final int totalFragments, final int uniqueFragments,
            final int sharedFragments, final int wildcardFragments) {
        ATypes = aTypes;
        BTypes = bTypes;
        CTypes = cTypes;
        TotalFragments = totalFragments;
        UniqueFragments = uniqueFragments;
        SharedFragments = sharedFragments;
        WildcardFragments = wildcardFragments;

        FittedFragments = uniqueFragments + sharedFragments + wildcardFragments;
        UnusedFragments = totalFragments - FittedFragments;
        PercentUnique = 1.0 * uniqueFragments / FittedFragments;
        PercentShared = 1.0 * sharedFragments / FittedFragments;
        PercentWildcard = 1.0 * wildcardFragments / FittedFragments;
    }

    public static CoverageQC create(
            int totalFragments, final HlaComplexCoverage winner)
    {
        List<HlaAllele> alleles = winner.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList());
        int aTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_A)).collect(Collectors.toSet()).size();
        int bTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_B)).collect(Collectors.toSet()).size();
        int cTypes = alleles.stream().filter(x -> x.Gene.equals(GENE_C)).collect(Collectors.toSet()).size();

        if (aTypes == 0 || bTypes == 0 || cTypes == 0)
        {
            LL_LOGGER.warn("  UNMATCHED_TYPE - {} A alleles, {} B alleles, {} C alleles", aTypes, bTypes, cTypes);
        }

        if (winner.WildCoverage > 0)
        {
            LL_LOGGER.warn("  WILDCARD_MATCH - winning solution contains wildcards");
        }

        return new CoverageQC(aTypes, bTypes, cTypes, totalFragments, winner.UniqueCoverage, winner.SharedCoverage, winner.WildCoverage);
    }

    public List<String> header()
    {
        return Lists.newArrayList(
                "aTypes", "bTypes", "cTypes", "unusedFragments", "fittedFragments",
                "percentUnique", "percentShared", "percentWildcard");
    }

    public List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(ATypes),
                String.valueOf(BTypes), String.valueOf(CTypes),
                String.valueOf(UnusedFragments),
                String.valueOf(FittedFragments),
                String.format("%.3f", PercentUnique),
                String.format("%.3f", PercentShared),
                String.format("%.3f", PercentWildcard));
    }

    public String toString()
    {
        return String.format("CoverageQC(aTypes=%d, bTypes=%d, cTypes=%d, unusedFragments=%d, uniqueFragments=%d, sharedFragments=%d, "
                        + "wildcardFragments=%d, fittedFragments=%d, percentUnique=%.3f, percentShared=%.3f, percentWildcard=%.3f)",
                ATypes, BTypes, CTypes, UnusedFragments, UniqueFragments, SharedFragments, FittedFragments,
                PercentUnique, PercentShared, PercentWildcard);
    }
}