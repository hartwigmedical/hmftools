package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;
import static java.lang.Math.round;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public final class HlaComplexCoverage implements Comparable<HlaComplexCoverage>
{
    public final int TotalCoverage;
    public final int UniqueCoverage;
    public final int SharedCoverage;
    public final int WildCoverage;
    
    private final List<HlaAlleleCoverage> mAlleleCoverage;

    public HlaComplexCoverage(int uniqueCoverage, int sharedCoverage, int wildCoverage, final List<HlaAlleleCoverage> alleleCoverage)
    {
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        mAlleleCoverage = alleleCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;
    }

    public List<HlaAlleleCoverage> getAlleleCoverage() { return mAlleleCoverage; }

    public final HlaComplexCoverage expandToSixAlleles()
    {
        return create(HlaAlleleCoverage.expand(mAlleleCoverage));
    }

    public final int homozygousAlleles()
    {
        int alleleCount = (int)mAlleleCoverage.stream().map(x -> x.Allele).distinct().count();
        return 6 - alleleCount;
    }

    @Override
    public int compareTo(final HlaComplexCoverage other)
    {
        if(TotalCoverage != other.TotalCoverage)
            return TotalCoverage > other.TotalCoverage ? 1 : -1;

        if(SharedCoverage != other.SharedCoverage)
            return SharedCoverage > other.SharedCoverage ? 1 : -1;

        if(WildCoverage != other.WildCoverage)
            return WildCoverage > other.WildCoverage ? 1 : -1;

        if(UniqueCoverage != other.UniqueCoverage)
            return UniqueCoverage > other.UniqueCoverage ? 1 : -1;

        return 0;
    }

    public String toString()
    {
        Set<HlaAllele> alleles = mAlleleCoverage.stream().map(x -> x.Allele).collect(Collectors.toSet());
        StringJoiner sj = new StringJoiner("\t");

        // CHECK what is joined?
        mAlleleCoverage.forEach(x -> sj.add(x.toString()));
        return String.format("%d\t%d\t%d\t%d\t%d\t%s",
                TotalCoverage, UniqueCoverage, SharedCoverage, WildCoverage, alleles.size(), sj.toString());
    }

    public static HlaComplexCoverage create(final List<HlaAlleleCoverage> alleles)
    {
        int unique = 0;
        double shared = 0.0;
        double wild = 0.0;
        for(HlaAlleleCoverage coverage : alleles)
        {
            unique += coverage.UniqueCoverage;
            shared += coverage.SharedCoverage;
            wild += coverage.WildCoverage;
        }

        final List<HlaAlleleCoverage> sortedAlleles = alleles.stream().collect(Collectors.toList());
        Collections.sort(sortedAlleles, new HlaAlleleCoverage.AlleleSorter());

        return new HlaComplexCoverage(unique, (int)round(shared), (int)round(wild), sortedAlleles);
    }

    public static String header()
    {
        return "totalCoverage\tuniqueCoverage\tsharedCoverage\twildCoverage\ttypes\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6";
    }

    public List<HlaAlleleCoverage> confirmUnique(final LilacConfig config)
    {
        List<HlaAlleleCoverage> unique = mAlleleCoverage.stream()
                .filter(x -> x.UniqueCoverage >= config.MinConfirmedUniqueCoverage).collect(Collectors.toList());

        Collections.sort(unique, Collections.reverseOrder());

        List<HlaAlleleCoverage> aList = takeN(unique.stream().filter(x -> x.Allele.Gene.equals("A")).collect(Collectors.toList()), 2);
        List<HlaAlleleCoverage> bList = takeN(unique.stream().filter(x -> x.Allele.Gene.equals("B")).collect(Collectors.toList()), 2);
        List<HlaAlleleCoverage> cList = takeN(unique.stream().filter(x -> x.Allele.Gene.equals("C")).collect(Collectors.toList()), 2);

        List<HlaAlleleCoverage> results = Lists.newArrayList();
        results.addAll(aList);
        results.addAll(bList);
        results.addAll(cList);
        Collections.sort(results, Collections.reverseOrder());
        return results;
    }

    private static List<HlaAlleleCoverage> takeN(final List<HlaAlleleCoverage> list, int n)
    {
        List<HlaAlleleCoverage> newList = Lists.newArrayList();

        for(int i = 0; i < min(list.size(), n); ++i)
        {
            newList.add(list.get(i));
        }

        return newList;
    }

    public final void writeToFile(final List<HlaComplexCoverage> coverage, final String fileName)
    {
        // TODO
        File file = new File(fileName);

        /*
        FilesKt.writeText$default((File) file, (String) (header() + "\n"), null, (int) 2, null);
        for(HlaComplexCoverage coverage : $receiver)
        {
            FilesKt.appendText$default((File) file, (String) (coverage.toString() + "\n"), null, (int) 2, null);
        }
        */
    }

    public static class TotalCoverageSorter implements Comparator<HlaComplexCoverage>
    {
        public int compare(final HlaComplexCoverage first, final HlaComplexCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage > second.TotalCoverage ? 1 : -1;

            return 0;
        }
    }

}
