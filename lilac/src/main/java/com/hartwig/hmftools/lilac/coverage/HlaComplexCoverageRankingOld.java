package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

/**
 * Complexes within [maxDistanceFromTopScore] fragments of the highest aligned fragment count are considered as possible solutions.
 * These are then ranked by:
 * - first prioritising solutions with the fewest alleles with wildcard matches,
 * - then choosing solutions with the most homozygous alleles,
 * - then choosing solution with the most common alleles,
 * - then choosing solution with the least recovered alleles
 * - finally choosing the solution with the lowest number.
 */

public class HlaComplexCoverageRankingOld implements Comparator<HlaComplexCoverage>
{
    private final int mMaxDistanceFromTopScore;
    private final List<HlaAllele> mCommon;
    private final List<HlaAllele> mRecovered;
    private final List<HlaAllele> mFavourites;

    public HlaComplexCoverageRankingOld(
            int maxDistanceFromTopScore, final List<HlaAllele> common, final List<HlaAllele> recovered, final List<HlaAllele> favourites)
    {
        mMaxDistanceFromTopScore = maxDistanceFromTopScore;
        mCommon = common;
        mRecovered = recovered;
        mFavourites = favourites;
    }

    public List<HlaComplexCoverage> candidateRanking(final List<HlaComplexCoverage> complexes)
    {
        int topScore = complexes.stream().mapToInt(x -> x.TotalCoverage).max().orElse(0);
        if (topScore == 0)
            return Lists.newArrayList();

        List<HlaComplexCoverage> results = complexes.stream()
                .filter(x -> x.TotalCoverage >= topScore - mMaxDistanceFromTopScore).collect(Collectors.toList());

        Collections.sort(results, this);
        return results;
    }

    public static final int homozygousAlleles(final HlaComplexCoverage coverage)
    {
        return EXPECTED_ALLELE_COUNT - coverage.homozygousCount();
    }

    public int compare(final HlaComplexCoverage o1, final HlaComplexCoverage o2)
    {
        int favouritesCount = favourites(o1) - favourites(o2);
        if (favouritesCount != 0)
            return -favouritesCount;

        int wildcardCount = wildcardCount(o1) - wildcardCount(o2);
        if (wildcardCount != 0)
            return wildcardCount;

        int homozygousCompare = homozygousAlleles(o1) - homozygousAlleles(o2);
        if (homozygousCompare != 0)
            return -homozygousCompare;

        int commonCountCompare = commonCount(o1) - commonCount(o2);
        if (commonCountCompare != 0)
            return -commonCountCompare;

        int recoveredCountCompare = recoveredCount(o1) - recoveredCount(o2);
        if (recoveredCountCompare != 0)
            return recoveredCountCompare;

        for (int i = 0; i < min(o1.getAlleleCoverage().size(), o2.getAlleleCoverage().size()); ++i)
        {
            HlaAllele o1Allele = o1.getAlleleCoverage().get(i).Allele;
            HlaAllele o2Allele = o2.getAlleleCoverage().get(i).Allele;

            int alleleCompare = o1Allele.compareTo(o2Allele);
            if (alleleCompare != 0)
                return alleleCompare;
        }

        return 0;
    }

    private int favourites(final HlaComplexCoverage coverage)
    {
        return (int)coverage.getAlleleCoverage().stream().filter(x -> mFavourites.contains(x.Allele)).count();
    }

    private static int wildcardCount(final HlaComplexCoverage coverage)
    {
        return (int)coverage.getAlleleCoverage().stream().filter(x -> x.WildCoverage > 0).count();
    }

    private final int commonCount(final HlaComplexCoverage coverage)
    {
        return (int)coverage.getAlleleCoverage().stream().filter(x -> mCommon.contains(x.Allele)).count();
    }

    private final int recoveredCount(final HlaComplexCoverage coverage)
    {
        return (int)coverage.getAlleleCoverage().stream().filter(x -> mRecovered.contains(x.Allele)).count();
    }
}
