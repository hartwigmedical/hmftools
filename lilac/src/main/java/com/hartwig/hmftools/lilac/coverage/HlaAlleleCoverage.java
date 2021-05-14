package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.hla.HlaAllele.contains;

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
        List<HlaAlleleCoverage> a = coverage.stream().filter(x -> x.Allele.Gene.equals("A")).collect(Collectors.toList());
        List<HlaAlleleCoverage> b = coverage.stream().filter(x -> x.Allele.Gene.equals("B")).collect(Collectors.toList());
        List<HlaAlleleCoverage> c = coverage.stream().filter(x -> x.Allele.Gene.equals("C")).collect(Collectors.toList());

        List<HlaAlleleCoverage> expandedCoverage = Lists.newArrayList();
        expandedCoverage.addAll(splitSingle(a));
        expandedCoverage.addAll(splitSingle(b));
        expandedCoverage.addAll(splitSingle(c));

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
        List<HlaAlleleCoverage> results = Lists.newArrayList();

        Map<HlaAllele,Integer> uniqueCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> combinedCoverageMap = Maps.newHashMap();
        Map<HlaAllele,Double> wildCoverageMap = Maps.newHashMap();

        for(FragmentAlleles fragment : fragmentSequences)
        {
            List<HlaAllele> fullAlleles = HlaAllele.dedup(fragment.getFull().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toList()));
            List<HlaAllele> partialAlleles = HlaAllele.dedup(fragment.getPartial().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toList()));
            List<HlaAllele> wildAlleles = HlaAllele.dedup(fragment.getWild().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toList()));

            if(fullAlleles.size() == 1 && partialAlleles.isEmpty())
            {
                increment(uniqueCoverageMap, fullAlleles.get(0), 1);
            }
            else
            {
                double contribution = 1.0 / (fullAlleles.size() + partialAlleles.size() + wildAlleles.size());
                fullAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                partialAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                wildAlleles.forEach(x -> increment(wildCoverageMap, x, contribution));
            }
        }

        List<HlaAllele> hlaAlleles = uniqueCoverageMap.keySet().stream().collect(Collectors.toList());
        combinedCoverageMap.keySet().stream().filter(x -> !contains(hlaAlleles, x)).forEach(x -> hlaAlleles.add(x));

        for (HlaAllele allele : hlaAlleles)
        {
            int uniqueCoverage = uniqueCoverageMap.entrySet().stream().filter(x -> x.getKey() == allele).mapToInt(x -> x.getValue()).findFirst().orElse(0);
            double combinedCoverage = combinedCoverageMap.entrySet().stream().filter(x -> x.getKey() == allele).mapToDouble(x -> x.getValue()).findFirst().orElse(0);
            double wildCoverage = wildCoverageMap.entrySet().stream().filter(x -> x.getKey() == allele).mapToDouble(x -> x.getValue()).findFirst().orElse(0);

            results.add(new HlaAlleleCoverage(allele, uniqueCoverage, combinedCoverage, wildCoverage));
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

            List<HlaAlleleCoverage> newCoverage = Lists.newArrayList(single, remainder);
            Collections.sort(newCoverage, new TotalCoverageSorter());
            return newCoverage;
        }

        return coverage;
    }

    public static List<HlaAllele> alleles(final List<HlaAlleleCoverage> coverage)
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
