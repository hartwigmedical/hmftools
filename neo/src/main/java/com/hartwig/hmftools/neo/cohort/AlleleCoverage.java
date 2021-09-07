package com.hartwig.hmftools.neo.cohort;

import java.util.List;

import com.google.common.collect.Lists;

public class AlleleCoverage
{
    public final String Allele;
    public final double CopyNumber;
    public final double VariantCount;

    public static final String DELIM = ",";
    public static final int EXPECTED_ALLELE_COUNT = 6;

    public AlleleCoverage(final String allele, double copyNumber, double variantCount)
    {
        Allele = allele;
        CopyNumber = copyNumber;
        VariantCount = variantCount;
    }

    public String gene() { return gene(Allele); }
    public static String gene(final String allele) { return allele.substring(0, 1); }

    public boolean isLost() { return CopyNumber < 0.5 || VariantCount > 0; }

    public String toString()
    {
        return String.format("allele(%s) copyNumber(%.2f) variants(%.1f)",
                Allele, CopyNumber, VariantCount);
    }

    public static AlleleCoverage fromCsv(final String data, int alleleIndex, int tumorCnIndex, final List<Integer> somaticVariantIndices)
    {
        // SomaticMissense>0|SomaticNonsenseOrFrameshift>0|SomaticSplice>0|SomaticInframeIndel>0
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("\\*", "").replaceAll(":", "");
        double variantCount = somaticVariantIndices.stream().mapToDouble(x -> Double.parseDouble(items[x])).sum();

        return new AlleleCoverage(allele, Double.parseDouble(items[tumorCnIndex]), variantCount);
    }

    public static List<Boolean> getGeneStatus(final List<AlleleCoverage> alleleCoverages)
    {
        final List<Boolean> geneLostStatus = Lists.newArrayListWithExpectedSize(EXPECTED_ALLELE_COUNT);

        for(AlleleCoverage alleleCoverage : alleleCoverages)
        {
            boolean geneLost = alleleCoverages.stream().filter(x -> x.gene().equals(alleleCoverage.gene())).anyMatch(x -> x.isLost());
            geneLostStatus.add(geneLost);
        }

        return geneLostStatus;
    }
}
