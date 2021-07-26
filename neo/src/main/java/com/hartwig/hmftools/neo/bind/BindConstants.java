package com.hartwig.hmftools.neo.bind;

import java.util.List;

import com.google.common.collect.Lists;

public final class BindConstants
{
    public static final List<Character> AMINO_ACIDS = Lists.newArrayList(
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

    public static final int MAX_PEPTIDE_LENGTH = 15;

    public static final double MIN_OBSERVED_AA_POS_FREQ = 0.001;

}
