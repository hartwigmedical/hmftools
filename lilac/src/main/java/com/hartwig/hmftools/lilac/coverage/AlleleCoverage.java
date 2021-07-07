package com.hartwig.hmftools.lilac.coverage;

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

public class AlleleCoverage implements Comparable<AlleleCoverage>
{
    public final HlaAllele Allele;
    public final int UniqueCoverage;
    public final double SharedCoverage;
    public final double WildCoverage;
    public final double TotalCoverage;

    public AlleleCoverage(final HlaAllele allele, int uniqueCoverage, double sharedCoverage, double wildCoverage)
    {
        Allele = allele;
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;
    }

    public int compareTo(final AlleleCoverage other)
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

    public static String toString(final List<AlleleCoverage> coverage)
    {
        StringJoiner sj = new StringJoiner(", ");
        coverage.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static List<AlleleCoverage> proteinCoverage(final List<FragmentAlleles> fragAlleles)
    {
        return buildCoverage(fragAlleles, false);
    }

    public static List<AlleleCoverage> groupCoverage(final List<FragmentAlleles> fragAlleles)
    {
        return buildCoverage(fragAlleles, true);
    }

    private static List<AlleleCoverage> buildCoverage(final List<FragmentAlleles> fragAlleles, boolean asAlleleGroup)
    {
        // attributes coverage counts to each allele amongst the set of alleles present across all fragments
        List<AlleleCoverage> results = Lists.newArrayList();

        Map<HlaAllele,Integer> uniqueCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> combinedCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> wildCoverageMap = Maps.newHashMap();

        // test whether a fragment has a single (unique) allele in its full list (with nothing in partial)
        // if not, split te contribution across the alleles into combined and wild

        for(FragmentAlleles fragment : fragAlleles)
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
            Integer uniqueCoverage = uniqueCoverageMap.get(allele);
            Double combinedCoverage = combinedCoverageMap.get(allele);
            Double wildCoverage = wildCoverageMap.get(allele);

            results.add(new AlleleCoverage(
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

    public static List<HlaAllele> coverageAlleles(final List<AlleleCoverage> coverage)
    {
        return coverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }

    public static class TotalCoverageSorter implements Comparator<AlleleCoverage>
    {
        // sorts by total coverage descending
        public int compare(final AlleleCoverage first, final AlleleCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage < second.TotalCoverage ? 1 : -1;

            return 0;
        }
    }

    public static class AlleleSorter implements Comparator<AlleleCoverage>
    {
        public int compare(final AlleleCoverage first, final AlleleCoverage second)
        {
            return first.Allele.toString().compareTo(second.Allele.toString());
        }
    }
}
