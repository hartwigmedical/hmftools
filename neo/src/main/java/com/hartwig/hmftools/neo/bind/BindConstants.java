package com.hartwig.hmftools.neo.bind;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public final class BindConstants
{
    public static final List<Character> AMINO_ACIDS = Lists.newArrayList(
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

    public static final int AMINO_ACID_COUNT = AMINO_ACIDS.size();

    public static final Map<Character,Integer> AMINO_ACID_INDICES = Maps.newHashMap();

    static
    {
        for(int i = 0; i < AMINO_ACIDS.size(); ++i)
        {
            AMINO_ACID_INDICES.put(AMINO_ACIDS.get(i), i);
        }
    }

    public static final int INVALID_AMINO_ACID = -1;

    public static int aminoAcidIndex(final char aminoAcid)
    {
        Integer index = AMINO_ACID_INDICES.get(aminoAcid);
        return index != null ? index : INVALID_AMINO_ACID;
    }

    public static final int REF_PEPTIDE_LENGTH = 12;
    public static final int ALLELE_POS_MAPPING_PEPTIDE_LENGTH = 9;
    public static final int REF_PEPTIDE_LEFT_FIXED_POS = 4;

    public static final double MIN_OBSERVED_AA_POS_FREQ = 0.005;

    public static final double AMINO_ACID_C_FREQ_ADJUST = 3;

    public static final double DEFAULT_MAX_AFFINITY = 20000;
    public static final double DEFAULT_BINDING_AFFINITY_LOW = 250;
    public static final double DEFAULT_BINDING_AFFINITY_HIGH = 500;

    public static final double DEFAULT_PEPTIDE_LENGTH_WEIGHT = 1000;
    public static final double DEFAULT_ALLELE_MOTIF_WEIGHT = 2000;
    public static final double DEFAULT_WEIGHT_EXPONENT = 1.5;

    public static final double DEFAULT_NOISE_PROB = 0.05;
    public static final double DEFAULT_NOISE_WEIGHT = 0.5;
}
