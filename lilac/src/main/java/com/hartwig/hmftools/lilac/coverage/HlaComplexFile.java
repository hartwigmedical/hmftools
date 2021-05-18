package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaComplexFile
{
    public static String header()
    {
        StringJoiner sb = new StringJoiner(DELIM);
        sb.add("score");
        sb.add("homozygousCount");
        sb.add("cohortFrequency");
        sb.add("totalCoverage");
        sb.add("uniqueCoverage");
        sb.add("sharedCoverage");
        sb.add("wildCoverage");
        sb.add("allele1");
        sb.add("allele2");
        sb.add("allele3");
        sb.add("allele4");
        sb.add("allele5");
        sb.add("allele6");
        return sb.toString();
    }

    public static void writeToFile(final List<HlaComplexCoverage> coverages, final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(header());
            writer.newLine();

            for(HlaComplexCoverage coverage : coverages)
            {
                writer.write(asString(coverage));
            }

            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
    }

    public static String asString(final HlaComplexCoverage coverage)
    {
        StringJoiner sj = new StringJoiner("\t");

        sj.add(String.format("%.2f", coverage.getScore()));
        sj.add(String.valueOf(coverage.homozygousCount()));
        sj.add(String.format("%.2f", coverage.cohortFrequencyTotal()));
        sj.add(String.valueOf(coverage.TotalCoverage));
        sj.add(String.valueOf(coverage.UniqueCoverage));
        sj.add(String.valueOf(coverage.SharedCoverage));
        sj.add(String.valueOf(coverage.WildCoverage));
        coverage.getAlleleCoverage().forEach(x -> sj.add(x.toString()));

        return sj.toString();
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

}
