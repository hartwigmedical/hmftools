package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;

import java.io.BufferedWriter;
import java.io.IOException;
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

    public List<HlaAlleleCoverage> confirmUnique(final LilacConfig config, final List<HlaAllele> confirmedGroupAlleles)
    {
        List<HlaAlleleCoverage> unique = mAlleleCoverage.stream()
                .filter(x -> x.UniqueCoverage >= config.MinConfirmedUniqueCoverage).collect(Collectors.toList());

        Collections.sort(unique, Collections.reverseOrder());

        List<HlaAlleleCoverage> results = Lists.newArrayList();

        // take at most 2 alleles for each gene, and at most 1 unique protein if more than 1 unique group is provided
        for(String gene : GENE_IDS)
        {
            List<HlaAlleleCoverage> geneCoverage = unique.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

            if(geneCoverage.isEmpty())
                continue;

            final List<HlaAllele> geneGroupAlleles = confirmedGroupAlleles.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList());

            int geneCount = 0;
            for(HlaAlleleCoverage coverage : geneCoverage)
            {
                if(geneGroupAlleles.size() > 1)
                {
                    // how many added already from this protein's group
                    int matchedGroupCount = (int)results.stream()
                            .filter(x -> geneGroupAlleles.contains(x.Allele.asAlleleGroup())
                                    && x.Allele.asAlleleGroup().equals(coverage.Allele.asAlleleGroup())).count();

                    if(matchedGroupCount >= 1)
                        continue;
                }

                results.add(coverage);

                ++geneCount;

                if(geneCount >= 2)
                    break;
            }
        }

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

    public static void writeToFile(final List<HlaComplexCoverage> coverage, final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(header());
            writer.newLine();

            writer.write(coverage.toString());
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
    }

    public static List<HlaAllele> parseCandidateCoverageData(final List<String> alleleDataList)
    {
        // convert allele coverage output back into the 6 candidate alleles
        // example: A*01:01[199,131,68,0]   A*02:01[162,100,62,0]   B*18:01[182,112,70,0]   B*38:01[165,92,73,0]    C*12:03[356,320,36,0]

        List<HlaAllele> rawAlleles = org.apache.commons.compress.utils.Lists.newArrayList();

        for(String alleleCoverageData : alleleDataList)
        {
            String alleleData = alleleCoverageData.replaceAll("\\[[0-9,]*]", "");
            rawAlleles.add(HlaAllele.fromString(alleleData));
        }

        List<HlaAllele> allAlleles = Lists.newArrayList();

        int aCount = 0;
        int bCount = 0;
        int cCount = 0;

        for(HlaAllele allele : rawAlleles)
        {
            if(allele.Gene.equals(GENE_A))
            {
                ++aCount;
                allAlleles.add(allele);
            }
            else if(allele.Gene.equals(GENE_B))
            {
                if(aCount == 1)
                {
                    ++aCount;
                    allAlleles.add(allAlleles.get(allAlleles.size() - 1));
                }

                ++bCount;
                allAlleles.add(allele);
            }
            else
            {
                if(bCount == 1)
                {
                    ++bCount;
                    allAlleles.add(allAlleles.get(allAlleles.size() - 1));
                }

                ++cCount;
                allAlleles.add(allele);
            }
        }

        if(cCount == 1)
            allAlleles.add(allAlleles.get(allAlleles.size() - 1));

        return allAlleles;
    }

    public static class TotalCoverageSorter implements Comparator<HlaComplexCoverage>
    {
        // sorts by total coverage descending
        public int compare(final HlaComplexCoverage first, final HlaComplexCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage < second.TotalCoverage ? 1 : -1;

            return 0;
        }
    }

}
