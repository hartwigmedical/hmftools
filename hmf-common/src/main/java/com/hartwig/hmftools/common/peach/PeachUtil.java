package com.hartwig.hmftools.common.peach;

import org.jetbrains.annotations.NotNull;

public class PeachUtil
{
    public static final String HAPLOTYPE_SEPARATOR = "_";
    public static final String DRUG_SEPARATOR = ";";

    public static final String HOMOZYGOUS_ZYGOSITY_STRING = "HOM";
    public static final String HETEROZYGOUS_ZYGOSITY_STRING = "HET";

    public static final String UNKNOWN_ALLELE_STRING = "Unresolved Haplotype";
    public static final String NOT_APPLICABLE_ZYGOSITY_STRING = "N/A";

    @NotNull
    public static String convertToHaplotypeString(@NotNull final String allele, int alleleCount)
    {
        if(allele.equals(UNKNOWN_ALLELE_STRING))
        {
            return allele;
        }
        else
        {
            return allele + HAPLOTYPE_SEPARATOR + convertToZygosityString(alleleCount);
        }
    }

    public static int convertZygosityToAlleleCount(@NotNull final String zygosity)
    {
        if(zygosity.equals(HETEROZYGOUS_ZYGOSITY_STRING))
        {
            return 1;
        }
        else if(zygosity.equals(HOMOZYGOUS_ZYGOSITY_STRING) || zygosity.equals(NOT_APPLICABLE_ZYGOSITY_STRING))
        {
            return 2;
        }
        else
        {
            throw new IllegalArgumentException(String.format("Could not convert zygosity %s to an allele count", zygosity));
        }
    }

    @NotNull
    public static String convertToZygosityString(int alleleCount)
    {
        if(alleleCount == 1)
        {
            return HETEROZYGOUS_ZYGOSITY_STRING;
        }
        else if (alleleCount == 2)
        {
            return HOMOZYGOUS_ZYGOSITY_STRING;
        }
        else
        {
            throw new IllegalArgumentException(String.format("Could not convert allele count %s to a zygosity", alleleCount));
        }
    }
}
