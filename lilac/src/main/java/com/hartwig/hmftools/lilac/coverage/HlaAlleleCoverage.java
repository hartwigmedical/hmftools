package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaAlleleCoverage implements Comparable<HlaAlleleCoverage>
{
    public final HlaAllele Allele;
    public final int UniqueCoverage;
    public final double SharedCoverage;
    public final double WildCoverage;
    public final double TotalCoverage;

    public HlaAlleleCoverage(final HlaAllele allele, int uniqueCoverage, double sharedCoverage, double wildCoverage)
    {
        Allele = allele;
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;
    }

    public int compareTo(final HlaAlleleCoverage other)
    {
        if(UniqueCoverage != other.UniqueCoverage)
            return UniqueCoverage > other.UniqueCoverage ? 1 : -1;

        if(SharedCoverage != other.SharedCoverage)
            return SharedCoverage > other.SharedCoverage ? 1 : -1;

        return 0;
    }

    public String toString()
    {
        return String.format("%s[%.0f,%d,%.0f,%.0f]",
                Allele.toString(), TotalCoverage, UniqueCoverage, SharedCoverage, WildCoverage);
    }

    public static String toString(final List<HlaAlleleCoverage> coverage)
    {
        StringJoiner sj = new StringJoiner(", ");
        coverage.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static List<HlaAlleleCoverage> expand(final List<HlaAlleleCoverage> coverage)
    {
        if(coverage.size() == EXPECTED_ALLELE_COUNT)
            return coverage;

        List<HlaAlleleCoverage> expandedCoverage = Lists.newArrayList();

        for(String gene : GENE_IDS)
        {
            List<HlaAlleleCoverage> geneCoverage = coverage.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

            if(geneCoverage.size() == 2)
                expandedCoverage.addAll(geneCoverage);
            else
                expandedCoverage.addAll(splitHomozygousCoverage(geneCoverage));
        }

        Collections.sort(expandedCoverage, new AlleleSorter());
        return expandedCoverage;
    }

    public static List<HlaAlleleCoverage> proteinCoverage(final List<FragmentAlleles> fragmentSequences)
    {
        return create(fragmentSequences, false);
    }

    public static List<HlaAlleleCoverage> groupCoverage(final List<FragmentAlleles> fragmentSequences)
    {
        return create(fragmentSequences, true);
    }

    public static List<HlaAlleleCoverage> create(final List<FragmentAlleles> fragmentSequences, boolean asAlleleGroup)
    {
        // attributes coverage counts to each allele amongst the set of alleles present across all fragments
        List<HlaAlleleCoverage> results = Lists.newArrayList();

        Map<HlaAllele,Integer> uniqueCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> combinedCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> wildCoverageMap = Maps.newHashMap();

        // test whether a fragment has a single (unique) allele in its full list (with nothing in partial)
        // if not, split te contribution across the alleles into combined and wild

        for(FragmentAlleles fragment : fragmentSequences)
        {
            Set<HlaAllele> fullAlleles = fragment.getFull().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toSet());
            Set<HlaAllele> wildAlleles = fragment.getWild().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toSet());

            if(fullAlleles.size() == 1 && wildAlleles.isEmpty())
            {
                increment(uniqueCoverageMap, fullAlleles.iterator().next(), 1);
            }
            else
            {
                double contribution = 1.0 / (fullAlleles.size() + wildAlleles.size());
                fullAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                wildAlleles.forEach(x -> increment(wildCoverageMap, x, contribution));
            }
        }

        // gather up the unique set of supported alleles, then tally their results into a single allele coverage container
        List<HlaAllele> hlaAlleles = uniqueCoverageMap.keySet().stream().collect(Collectors.toList());
        combinedCoverageMap.keySet().stream().filter(x -> !hlaAlleles.contains(x)).forEach(x -> hlaAlleles.add(x));

        for (HlaAllele allele : hlaAlleles)
        {
            // int uniqueCoverage = uniqueCoverageMap.entrySet().stream().filter(x -> x.getKey() == allele).mapToInt(x -> x.getValue()).findFirst().orElse(0);
            Integer uniqueCoverage = uniqueCoverageMap.get(allele);
            Double combinedCoverage = combinedCoverageMap.get(allele);
            Double wildCoverage = wildCoverageMap.get(allele);

            results.add(new HlaAlleleCoverage(
                    allele,
                    uniqueCoverage != null ? uniqueCoverage : 0,
                    combinedCoverage != null ? combinedCoverage : 0,
                    wildCoverage != null ? wildCoverage : 0));
        }

        Collections.sort(results, Collections.reverseOrder());
        return results;
    }

    private static void increment(final Map<HlaAllele,Integer> map, final HlaAllele allele, int value)
    {
        Map.Entry<HlaAllele,Integer> entry = map.entrySet().stream().filter(x -> x.getKey() == allele).findFirst().orElse(null);

        if(entry == null)
        {
            map.put(allele, value);
            return;
        }

        entry.setValue(entry.getValue() + value);
    }

    private static void increment(final Map<HlaAllele,Double> map, final HlaAllele allele, double value)
    {
        Map.Entry<HlaAllele,Double> entry = map.entrySet().stream().filter(x -> x.getKey() == allele).findFirst().orElse(null);

        if(entry == null)
        {
            map.put(allele, value);
            return;
        }

        entry.setValue(entry.getValue() + value);
    }

    private static List<HlaAlleleCoverage> splitHomozygousCoverage(final List<HlaAlleleCoverage> coverage)
    {
        if (coverage.size() != 1)
            return coverage;

        HlaAlleleCoverage single = coverage.get(0);
        HlaAlleleCoverage first = new HlaAlleleCoverage(
                single.Allele,
                single.UniqueCoverage / 2,
                single.SharedCoverage / 2,
                single.WildCoverage / 2);

        HlaAlleleCoverage remainder = new HlaAlleleCoverage(single.Allele,
                single.UniqueCoverage - first.UniqueCoverage,
                single.SharedCoverage - first.SharedCoverage,
                single.WildCoverage - single.WildCoverage);

        List<HlaAlleleCoverage> newCoverage = Lists.newArrayList(first, remainder);
        return newCoverage;
    }

    public static List<HlaAllele> coverageAlleles(final List<HlaAlleleCoverage> coverage)
    {
        return coverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }

    public static class TotalCoverageSorter implements Comparator<HlaAlleleCoverage>
    {
        // sorts by total coverage descending
        public int compare(final HlaAlleleCoverage first, final HlaAlleleCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage < second.TotalCoverage ? 1 : -1;

            return 0;
        }
    }

    public static class AlleleSorter implements Comparator<HlaAlleleCoverage>
    {
        public int compare(final HlaAlleleCoverage first, final HlaAlleleCoverage second)
        {
            return first.Allele.toString().compareTo(second.Allele.toString());
        }
    }
}
