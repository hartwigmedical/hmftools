package com.hartwig.hmftools.finding;

public class ChromosomeUtil
{
    public static String normalize(String chromosome) {
        return chromosome.replace("chr", "");
    }
}
