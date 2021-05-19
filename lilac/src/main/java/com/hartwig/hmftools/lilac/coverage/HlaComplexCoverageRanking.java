package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.abs;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.FREQUENCY_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.HOMOZYGOUS_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.RECOVERY_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.TOTAL_COVERAGE_DENOM;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.cohort.CohortFrequency;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

/* Rank candidates by:
    - first calculating bonuses and penalties for being homozygous and common in the provided cohort
    - ensure known frameshift alleles have been given unique fragment counts for each fragment supporting the indel
    - calculate a score as TotalCoverage + Adjusted Homozygous Count + Adjusted Cohort-Frequency Count
    - order by score
    - where 2 or more solutions have the same score, take the numerically (alphabetically) lowest solution

    Prevous ranking scheme:
    1. any solution with a 'favourite', which means the known frameshift allele
    2. fewest wildcard matches
    3. most homozygous alleles
    4. most common alleles
    5. least recovered alleles (seems redundant but I guess there could be a solution with common alleles, and another with the same number of common alleles, but some of them needing recovery?)
    6. lowest numbered alleles (sorted alphabetically)
* */

public class HlaComplexCoverageRanking
{
    private final int mMaxScoreDifference;
    private final ReferenceData mRefData;

    public HlaComplexCoverageRanking(final int maxScoreDifference, final ReferenceData refData)
    {
        mMaxScoreDifference = maxScoreDifference;
        mRefData = refData;
    }

    public List<HlaComplexCoverage> rankCandidates(final List<HlaComplexCoverage> complexes, final List<HlaAllele> recoveredAlleles)
    {
        if(complexes.isEmpty())
            return complexes;

        for(HlaComplexCoverage complexCoverage : complexes)
        {
            calcCohortFrequency(complexCoverage);
            calcRecoveryPenalty(complexCoverage, recoveredAlleles);
            calcComplexScore(complexCoverage);
        }

        if(mMaxScoreDifference == 0)
        {
            Collections.sort(complexes, new ComplexCoverageSorter());
            return complexes;
        }

        double topScore = complexes.stream().mapToDouble(x -> x.getScore()).max().orElse(0);
        int topCoverage = complexes.stream().mapToInt(x -> x.TotalCoverage).max().orElse(0);
        double inclusionThreshold = min(topScore - mMaxScoreDifference, topScore - mMaxScoreDifference * topCoverage / TOTAL_COVERAGE_DENOM);

        List<HlaComplexCoverage> results = complexes.stream()
                .filter(x -> x.getScore() >= inclusionThreshold).collect(Collectors.toList());

        Collections.sort(results, new ComplexCoverageSorter());
        return results;
    }

    private void calcRecoveryPenalty(final HlaComplexCoverage complexCoverage, final List<HlaAllele> recoveredAlleles)
    {
        int recoveredCount = (int)complexCoverage.getAlleles().stream().filter(x -> recoveredAlleles.contains(x)).count();
        complexCoverage.setRecoveredCount(recoveredCount);
    }

    private void calcCohortFrequency(final HlaComplexCoverage complexCoverage)
    {
        final CohortFrequency cohortFrequency = mRefData.getAlleleFrequencies();

        double cohortFrequencyTotal = 0;
        final List<HlaAllele> alleles = complexCoverage.getAlleles();

        for(HlaAllele allele : alleles)
        {
            double frequency = cohortFrequency.getAlleleFrequency(allele);
            cohortFrequencyTotal += log10(max(frequency, 0.0001));
        }

        complexCoverage.setCohortFrequencyTotal(cohortFrequencyTotal);
    }

    private void calcComplexScore(final HlaComplexCoverage complexCoverage)
    {
        int totalCoverage = complexCoverage.TotalCoverage;
        double adjustedCoverageFactor = totalCoverage / (double)TOTAL_COVERAGE_DENOM;

        double score = totalCoverage
                + complexCoverage.cohortFrequencyTotal() * FREQUENCY_SCORE_PENALTY * adjustedCoverageFactor
                + complexCoverage.homozygousCount() * HOMOZYGOUS_SCORE_PENALTY * adjustedCoverageFactor
                + complexCoverage.recoveredCount() * RECOVERY_SCORE_PENALTY * adjustedCoverageFactor;

        complexCoverage.setScore(score);
    }

    public static class ComplexCoverageSorter implements Comparator<HlaComplexCoverage>
    {
        // sorts by score then numerically if required
        public int compare(final HlaComplexCoverage first, final HlaComplexCoverage second)
        {
            if(abs(first.getScore() - second.getScore()) > 0.0001)
                return first.getScore() < second.getScore() ? 1 : -1;

            final List<HlaAllele> alleles1 = first.getAlleles();
            final List<HlaAllele> alleles2 = second.getAlleles();
            for (int i = 0; i < min(alleles1.size(), alleles2.size()); ++i)
            {
                HlaAllele o1Allele = alleles1.get(i);
                HlaAllele o2Allele = alleles2.get(i);

                int alleleCompare = o1Allele.compareTo(o2Allele);
                if (alleleCompare != 0)
                    return alleleCompare;
            }

            return 0;
        }
    }


}
