package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.abs;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.FREQUENCY_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.HOMOZYGOUS_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_POPULATION_FREQUENCY;
import static com.hartwig.hmftools.lilac.LilacConstants.RECOVERY_SCORE_PENALTY;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.cohort.CohortFrequency;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import com.google.common.collect.Lists;

/* Rank candidates by:
    - first calculating bonuses and penalties for being homozygous and common in the provided cohort
    - ensure known frameshift alleles have been given unique fragment counts for each fragment supporting the indel
    - calculate a score as TotalCoverage + Adjusted Homozygous Count + Adjusted Cohort-Frequency Count
    - order by score
    - where 2 or more solutions have the same score, take the numerically (alphabetically) lowest solution
*/

public class ComplexCoverageRanking
{
    private final double mMaxScoreDifference;
    private final ReferenceData mRefData;

    public ComplexCoverageRanking(double maxScoreDifference, final ReferenceData refData)
    {
        mMaxScoreDifference = maxScoreDifference;
        mRefData = refData;
    }

    public List<ComplexCoverage> rankCandidates(
            final List<ComplexCoverage> complexes, final List<HlaAllele> recoveredAlleles, final List<HlaSequenceLoci> sequences)
    {
        if(complexes.isEmpty())
            return complexes;

        for(ComplexCoverage complexCoverage : complexes)
        {
            calcScore(complexCoverage, recoveredAlleles, sequences);
        }

        if(mMaxScoreDifference == 0)
        {
            Collections.sort(complexes, new ComplexCoverageSorter());
            return complexes;
        }

        double topScore = complexes.stream().mapToDouble(x -> x.getScore()).max().orElse(0);
        int topCoverage = complexes.stream().mapToInt(x -> x.TotalCoverage).max().orElse(0);
        double inclusionThreshold = topScore - mMaxScoreDifference * topCoverage;

        // initial cull before sort
        double baseThreshold = topScore - 0.25 * topCoverage;
        List<ComplexCoverage> candidateResults = complexes.stream().filter(x -> x.getScore() >= baseThreshold).collect(Collectors.toList());
        Collections.sort(candidateResults, new ComplexCoverageSorter());

        List<ComplexCoverage> results = Lists.newArrayList();

        for(ComplexCoverage complexCoverage : candidateResults)
        {
            if(complexCoverage.getScore() >= inclusionThreshold || results.size() < 2)
            {
                results.add(complexCoverage);
            }
            else
            {
                break;
            }
        }

        return results;
    }

    private void calcScore(final ComplexCoverage complexCoverage, final List<HlaAllele> recoveredAlleles, final List<HlaSequenceLoci> sequences)
    {
        calcCohortFrequency(complexCoverage);
        calcRecoveryPenalty(complexCoverage, recoveredAlleles);
        calcWildcardPenalty(complexCoverage, sequences);
        calcComplexScore(complexCoverage);
    }

    private void calcRecoveryPenalty(final ComplexCoverage complexCoverage, final List<HlaAllele> recoveredAlleles)
    {
        int recoveredCount = (int)complexCoverage.getAlleles().stream()
                .filter(x -> recoveredAlleles.contains(x))
                .count();
        complexCoverage.setRecoveredCount(recoveredCount);
    }

    private void calcWildcardPenalty(final ComplexCoverage complexCoverage, final List<HlaSequenceLoci> sequences)
    {
        if(sequences.isEmpty())
            return;

        int wildcardCount = 0;

        for(HlaAllele allele : complexCoverage.getAlleles())
        {
            if(!allele.hasWildcards())
                continue;

            HlaSequenceLoci sequenceLoci = sequences.stream().filter(x -> x.Allele.equals(allele)).findFirst().orElse(null);

            if(sequenceLoci != null)
                wildcardCount += sequenceLoci.wildcardCount();
        }

        complexCoverage.setWildcardCount(wildcardCount);
    }

    private void calcCohortFrequency(final ComplexCoverage complexCoverage)
    {
        final CohortFrequency cohortFrequency = mRefData.getAlleleFrequencies();

        double cohortFrequencyTotal = 0;
        final List<HlaAllele> alleles = complexCoverage.getAlleles();

        for(HlaAllele allele : alleles)
        {
            double frequency = cohortFrequency.getAlleleFrequency(allele);
            double cohortPenalty = log10(max(frequency, MIN_POPULATION_FREQUENCY));

            cohortFrequencyTotal += cohortPenalty;

            if(complexCoverage.isHomozygous(allele))
                cohortFrequencyTotal += cohortPenalty;
        }

        complexCoverage.setCohortFrequencyTotal(cohortFrequencyTotal);
    }

    private void calcComplexScore(final ComplexCoverage complexCoverage)
    {
        int totalCoverage = complexCoverage.TotalCoverage;

        double score = totalCoverage
                + complexCoverage.cohortFrequencyTotal() * FREQUENCY_SCORE_PENALTY * totalCoverage
                + complexCoverage.homozygousCount() * HOMOZYGOUS_SCORE_PENALTY * totalCoverage
                - complexCoverage.recoveredCount() * RECOVERY_SCORE_PENALTY * totalCoverage;

        complexCoverage.setScore(score);
    }

    public static class ComplexCoverageSorter implements Comparator<ComplexCoverage>
    {
        // sorts by score then numerically if required
        public int compare(final ComplexCoverage first, final ComplexCoverage second)
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
