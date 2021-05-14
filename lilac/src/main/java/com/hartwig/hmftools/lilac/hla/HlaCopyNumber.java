package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

public class HlaCopyNumber
{
    public final HlaAllele Allele;
    public final double CopyNumber;

    public HlaCopyNumber(final HlaAllele allele, double copyNumber)
    {
        Allele = allele;
        CopyNumber = copyNumber;
    }

    public static List<HlaCopyNumber> alleleCopyNumber(final List<HlaAllele> winners)
    {
        return winners.stream().map(x -> new HlaCopyNumber(x, 0)).collect(Collectors.toList());
    }

    public static List<HlaCopyNumber> alleleCopyNumber(
            final List<HlaAllele> winners, final String geneCopyNumberFile, final HlaComplexCoverage tumorCoverage)
    {
        if(geneCopyNumberFile.isEmpty() || tumorCoverage.getAlleleCoverage().isEmpty())
            return winners.stream().map(x -> new HlaCopyNumber(x, 0.0)).collect(Collectors.toList());

        try
        {
            List<GeneCopyNumber> geneCopyNumberList = GeneCopyNumberFile.read(geneCopyNumberFile);

            List<HlaAlleleCoverage> alleleCoverage = tumorCoverage.getAlleleCoverage();

            List<GeneCopyNumber> hlaGeneCopyNumbers = geneCopyNumberList.stream()
                    .filter(x -> HLA_GENES.contains(x.gene())).collect(Collectors.toList());

            List<HlaCopyNumber> hlaCopyNumbers = Lists.newArrayList();

            for(String gene : HLA_GENES)
            {
                String geneId = gene.substring(gene.length() - 1);
                hlaCopyNumbers.addAll(alleleCopyNumber(
                        geneCopyNumberList.stream().filter(x -> x.gene().equals(gene)).findFirst().orElse(null),
                        alleleCoverage.stream().filter(x -> x.Allele.Gene.equals(geneId)).collect(Collectors.toList())));
            }

            return hlaCopyNumbers;
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read gene copy number file({}): {}", geneCopyNumberFile, e.toString());
            return winners.stream().map(x -> new HlaCopyNumber(x, 0.0)).collect(Collectors.toList());
        }
    }

    private static List<HlaCopyNumber> alleleCopyNumber(final GeneCopyNumber geneCopyNumber, final List<HlaAlleleCoverage> alleleCoverage)
    {
        if(geneCopyNumber == null || alleleCoverage.size() != 2)
            return Lists.newArrayList();

        double minor = geneCopyNumber.minMinorAlleleCopyNumber();
        double major = geneCopyNumber.minCopyNumber() - minor;

        HlaAllele majorAllele = alleleCoverage.get(0).TotalCoverage >= alleleCoverage.get(1).TotalCoverage ?
                alleleCoverage.get(0).Allele : alleleCoverage.get(1).Allele;

        HlaAllele minorAllele = alleleCoverage.get(0).Allele == majorAllele ?
                alleleCoverage.get(1).Allele : alleleCoverage.get(0).Allele;

        return Lists.newArrayList(new HlaCopyNumber(majorAllele, major), new HlaCopyNumber(minorAllele, minor));
    }

}
