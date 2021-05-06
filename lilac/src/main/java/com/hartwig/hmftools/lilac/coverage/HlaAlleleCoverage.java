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
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

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
        // CHECK - really store genes as just A, B and C?
        List<HlaAlleleCoverage> a = coverage.stream().filter(x -> x.Allele.Gene.equals("A")).collect(Collectors.toList());
        List<HlaAlleleCoverage> b = coverage.stream().filter(x -> x.Allele.Gene.equals("B")).collect(Collectors.toList());
        List<HlaAlleleCoverage> c = coverage.stream().filter(x -> x.Allele.Gene.equals("B")).collect(Collectors.toList());

        List<HlaAlleleCoverage> expandedCoverage = Lists.newArrayList();
        expandedCoverage.addAll(splitSingle(a));
        expandedCoverage.addAll(splitSingle(b));
        expandedCoverage.addAll(splitSingle(c));

        // CHECK
        Collections.sort(expandedCoverage, new AlleleSorter());
        return expandedCoverage;
    }

    public static List<HlaAlleleCoverage> proteinCoverage(final List<FragmentAlleles> fragmentSequences)
    {
        return create(fragmentSequences);
    }

    public static List<HlaAlleleCoverage> groupCoverage(final List<FragmentAlleles> fragmentSequences)
    {
        // TODO - check impact from orig
        // return create(fragmentSequences) { it.asAlleleGroup() }
        // return create(fragmentSequences).stream().map(x -> x.) { it.asAlleleGroup() }
        return Lists.newArrayList();
    }

    public static List<HlaAlleleCoverage> create(final List<FragmentAlleles> fragmentSequences) // pass in functor? final Function1<? super HlaAllele, HlaAllele> type
    {
        List<HlaAlleleCoverage> results = Lists.newArrayList();

        Map<HlaAllele,Integer> uniqueCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> combinedCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> wildCoverageMap = Maps.newHashMap();

        // Counts
        for(FragmentAlleles fragment : fragmentSequences)
        {
            Set<HlaAllele> fullAlleles = fragment.getFull().stream().collect(Collectors.toSet());
            Set<HlaAllele> partialAlleles = fragment.getPartial().stream().collect(Collectors.toSet());
            Set<HlaAllele> wildAlleles = fragment.getWild().stream().collect(Collectors.toSet());

            if(fullAlleles.size() == 1 && partialAlleles.isEmpty())
            {
                increment(uniqueCoverageMap, fullAlleles.iterator ().next(), 1);
            }
            else
            {
                double contribution = 1.0 / (fullAlleles.size() + partialAlleles.size() + wildAlleles.size());
                fullAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                partialAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                wildAlleles.forEach(x -> increment(wildCoverageMap, x, contribution));
            }
        }

        // CHECK unique alleles are added
        Set<HlaAllele> hlaAlleles = Sets.newHashSet();
        hlaAlleles.addAll(uniqueCoverageMap.keySet());
        hlaAlleles.addAll(combinedCoverageMap.keySet());

        for (HlaAllele allele : hlaAlleles)
        {
            int uniqueCoverage = uniqueCoverageMap.containsKey(allele) ? uniqueCoverageMap.get(allele) : 0;
            double combinedCoverage = combinedCoverageMap.containsKey(allele) ? combinedCoverageMap.get(allele) : 0;
            double wildCoverage = wildCoverageMap.containsKey(allele) ? wildCoverageMap.get(allele) : 0;

            results.add(new HlaAlleleCoverage(allele, uniqueCoverage, combinedCoverage, wildCoverage));
        }

        Collections.sort(results, Collections.reverseOrder());
        return results;
    }

    private static void increment(final Map<HlaAllele,Integer> map, final HlaAllele allele, int value)
    {
        Integer count = map.get(allele);
        map.put(allele, count != null ? count + value : value);
    }

    private static void increment(final Map<HlaAllele,Double> map, final HlaAllele allele, double value)
    {
        Double count = map.get(allele);
        map.put(allele, count != null ? count + value : value);
    }

    private static List<HlaAlleleCoverage> splitSingle(final List<HlaAlleleCoverage> coverage)
    {
        if (coverage.size() == 1)
        {
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

            // CHECK ordering
            List<HlaAlleleCoverage> newCoverage = Lists.newArrayList(single, remainder);
            Collections.sort(newCoverage, new TotalCoverageSorter());
            return newCoverage;

            /*
            if(single.TotalCoverage > remainder.TotalCoverage)
                return Lists.newArrayList(single, remainder);
            else
                return Lists.newArrayList(remainder, single);
             */

            // return Lists.newArrayList(first, remainder).sortedBy { it.totalCoverage }.reversed()
        }

        return coverage;
    }

    public static List<HlaAllele> alleles(final List<HlaAlleleCoverage> coverage)
    {
        return coverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }

    private static class TotalCoverageSorter implements Comparator<HlaAlleleCoverage>
    {
        public int compare(final HlaAlleleCoverage first, final HlaAlleleCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage > second.TotalCoverage ? 1 : -1;

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
