package com.hartwig.hmftools.neo.bind;

import java.util.List;

import com.google.common.collect.Lists;

public final class BindCommon
{
    public static final String FLD_ALLELE = "Allele";
    public static final String FLD_PEPTIDE = "Peptide";
    public static final String FLD_POSITION = "Position";
    public static final String FLD_SOURCE = "Source";
    public static final String FLD_UP_FLANK = "UpFlank";
    public static final String FLD_DOWN_FLANK = "DownFlank";
    public static final String FLD_AMINO_ACID = "AminoAcid";
    public static final String FLD_PEPTIDE_LEN = "PeptideLength";

    public static final String FLD_AFFINITY = "Affinity";
    public static final String FLD_PRED_AFFINITY = "PredictedAffinity";
    public static final String FLD_PRES_SCORE = "PresentationScore";

    public static final String DATA_TYPE_POS_WEIGHTS = "PosWeights";
    public static final String DATA_TYPE_BIND_COUNTS = "BindCounts";

    public static final List<String> COUNT_DATA_TYPES = Lists.newArrayList(
            DATA_TYPE_BIND_COUNTS, "NoiseCounts", "PeptideLengthWeighted", "AlleleMotifWeighted");

    public static final String DELIM = ",";
    public static final String RANDOM_SOURCE = "Random";
}
