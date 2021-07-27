package com.hartwig.hmftools.neo.bind;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public final class BindConstants
{
    public static final List<Character> AMINO_ACIDS = Lists.newArrayList(
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

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

    public static final int MAX_PEPTIDE_LENGTH = 15;

    public static final double MIN_OBSERVED_AA_POS_FREQ = 0.001;

}
