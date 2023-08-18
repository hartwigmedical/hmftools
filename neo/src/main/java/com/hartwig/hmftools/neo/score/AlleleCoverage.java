package com.hartwig.hmftools.neo.score;

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
}
