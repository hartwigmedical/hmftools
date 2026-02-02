package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.file.FileLock;
import com.hartwig.hmftools.lilac.GeneSelector;
import com.hartwig.hmftools.lilac.hla.HlaGene;

import org.jetbrains.annotations.Nullable;

public final class HlaComplexFile
{
    private HlaComplexFile() {}

    public static String header()
    {
        StringJoiner sb = new StringJoiner(TSV_DELIM);
        sb.add("Genes");
        sb.add("Score");
        sb.add("ComplexityPenalty");
        sb.add("Complexity");
        sb.add("HomozygousCount");
        sb.add("CohortFrequencyPenalty");
        sb.add("CohortFrequency");
        sb.add("RecoveryCount");
        sb.add("WildcardCount");
        sb.add("TotalCoverage");
        sb.add("UniqueCoverage");
        sb.add("SharedCoverage");
        sb.add("WildCoverage");
        sb.add("AlleleInfo");

        return sb.toString();
    }

    public static String infoHeader(final GeneSelector genes)
    {
        StringJoiner sb = new StringJoiner(TSV_DELIM);
        sb.add("Score");
        sb.add("ComplexityPenalty");
        sb.add("Complexity");
        sb.add("HomozygousCount");
        sb.add("CohortFrequencyPenalty");
        sb.add("CohortFrequency");
        sb.add("RecoveryCount");
        sb.add("WildcardCount");
        sb.add("TotalCoverage");
        sb.add("UniqueCoverage");
        sb.add("SharedCoverage");
        sb.add("WildCoverage");
        for(HlaGene gene : genes.genes())
        {
            sb.add(gene.shortName() + "_1");
            sb.add(gene.shortName() + "_2");
        }

        return sb.toString();
    }

    public static void writeToFile(final String fileName, final GeneSelector currentGenes, final Iterable<ComplexCoverage> coverages)
    {
        File file = new File(fileName);
        List<String> existingRecords = Lists.newArrayList();
        try(FileLock fileLock = FileLock.create(file))
        {
            BufferedReader reader = fileLock.getBufferedReader();
            reader.readLine();
            String line = reader.readLine();
            while(line != null)
            {
                String genesField = line.split(TSV_DELIM)[0];
                GeneSelector genes;
                try
                {
                    genes = GeneSelector.valueOf(genesField);
                }
                catch(IllegalArgumentException e)
                {
                    genes = null;
                }

                if(genes != null && genes != currentGenes)
                    existingRecords.add(line);

                line = reader.readLine();
            }

            fileLock.clear();
            BufferedWriter writer = fileLock.getBufferedWriter();
            writer.write(header());
            writer.newLine();
            for(String existingRecord : existingRecords)
            {
                writer.write(existingRecord);
                writer.newLine();
            }

            for(ComplexCoverage coverage : coverages)
            {
                writer.write(asString(currentGenes, coverage));
                writer.newLine();
            }

            writer.flush();
        }
        catch(Exception e)
        {
            LL_LOGGER.error("failed to update {}: {}", fileName, e.toString());
        }
    }

    public static String asString(@Nullable final GeneSelector genes, final ComplexCoverage coverage)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(genes != null)
            sj.add(genes.name());

        sj.add(String.format("%.2f", coverage.getScore()));
        sj.add(String.format("%.2f", coverage.getComplexityPenalty()));
        sj.add(String.format("%.2f", coverage.getComplexity()));
        sj.add(String.valueOf(coverage.homozygousCount()));
        sj.add(String.format("%.2f", coverage.getCohortFrequencyPenalty()));
        sj.add(String.format("%.2f", coverage.cohortFrequencyTotal()));
        sj.add(String.valueOf(coverage.recoveredCount()));
        sj.add(String.valueOf(coverage.wildcardCount()));
        sj.add(String.valueOf(coverage.TotalCoverage));
        sj.add(String.valueOf(coverage.UniqueCoverage));
        sj.add(String.valueOf(coverage.SharedCoverage));
        sj.add(String.valueOf(coverage.WildCoverage));
        if(genes != null)
        {
            StringJoiner alleleInfoBuilder = new StringJoiner(";");
            coverage.getAlleleCoverage().forEach(x -> alleleInfoBuilder.add(x.toString()));
            sj.add(alleleInfoBuilder.toString());
        }
        else
        {
            coverage.getAlleleCoverage().forEach(x -> sj.add(x.toString()));
        }

        return sj.toString();
    }

    public static void writeFragmentAssignment(
            final BufferedWriter writer, final Iterable<ComplexCoverage> coverages, final List<FragmentAlleles> fragAlleles)
            throws IOException
    {
        StringJoiner sb = new StringJoiner(TSV_DELIM);
        sb.add("Complex");
        sb.add("FragmentId");
        sb.add("FragmentCoords");
        sb.add("Type");
        sb.add("FullAlleles");
        sb.add("WildAlleles");
        writer.write(sb.toString());
        writer.newLine();

        for(ComplexCoverage complexCoverage : coverages)
        {
            StringJoiner complexAlleles = new StringJoiner(ITEM_DELIM);
            complexCoverage.getAlleles().forEach(x -> complexAlleles.add(x.toString()));
            String complexStr = complexAlleles.toString();

            List<FragmentAlleles> complexFrags = FragmentAlleles.filter(fragAlleles, complexCoverage.getAlleles());

            for(FragmentAlleles fragAllele : complexFrags)
            {
                StringJoiner fullAlleles = new StringJoiner(ITEM_DELIM);
                fragAllele.getFull().forEach(x -> fullAlleles.add(x.toString()));

                StringJoiner wildAlleles = new StringJoiner(ITEM_DELIM);
                fragAllele.getWild().forEach(x -> wildAlleles.add(x.toString()));

                String type = (fragAllele.getFull().size() == 1 && fragAllele.getWild().isEmpty()) ? "UNIQUE_FULL" :
                        (fragAllele.getFull().isEmpty() && fragAllele.getWild().size() == 1) ? "UNIQUE_WILD" : "SHARED";

                writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s",
                        complexStr, fragAllele.getFragment().id(), fragAllele.getFragment().readInfo(),
                        type, fullAlleles, wildAlleles));

                writer.newLine();
            }
        }
    }

    /* unused and Class-I specific
    public static List<HlaAllele> parseCandidateCoverageData(final List<String> alleleDataList)
    {
        // convert allele coverage output back into the 6 candidate alleles
        // example: A*01:01[199,131,68,0]   A*02:01[162,100,62,0]   B*18:01[182,112,70,0]   B*38:01[165,92,73,0]    C*12:03[356,320,36,0]

        List<HlaAllele> rawAlleles = new ArrayList<>();

        for(String alleleCoverageData : alleleDataList)
        {
            String alleleData = alleleCoverageData.replaceAll("\\[[0-9,]*]", "");
            rawAlleles.add(HlaAllele.fromString(alleleData));
        }

        List<HlaAllele> allAlleles = new ArrayList<>();

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
    */
}
