package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.redux.JitterModelParams.REPEAT_UNIT_3_PLUS_LABEL;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.JitterModelParams;

public final class JitterConstants
{
    // HP = homopolymer, DN = dinucleotide, 3P = 3+ length repeats

    // Illumina defaults
    public static final JitterVariableParams ILLUMINA_WGS_HP = new JitterVariableParams(
            0.105, 0.115, 0.125, 0.04);

    public static final JitterVariableParams ILLUMINA_WGS_DN = new JitterVariableParams(
            0.11, 0.15, 0.19, 0.04);

    public static final JitterVariableParams ILLUMINA_WGS_3P = new JitterVariableParams(
            0.14, 0.18, 0.22, 0.02);

    // Illumina high-depth defaults
    public static double DEFAULT_HD_JITTER_SCALE_GRADIENT = 0.06;

    public static final JitterVariableParams ILLUMINA_HIGH_DEPTH_HP = new JitterVariableParams(
            0.1, 0.15, 0.2, DEFAULT_HD_JITTER_SCALE_GRADIENT);

    public static final JitterVariableParams ILLUMINA_HIGH_DEPTH_DN = new JitterVariableParams(
            0.15, 0.2, 0.25, DEFAULT_HD_JITTER_SCALE_GRADIENT);

    public static final JitterVariableParams ILLUMINA_HIGH_DEPTH_3P = new JitterVariableParams(
            0.2, 0.25, 0.3, DEFAULT_HD_JITTER_SCALE_GRADIENT);

    // Ultima defaults
    public static final JitterVariableParams ULTIMA_WGS_HP = new JitterVariableParams(
            0.26, 0.29, 0.32, 0.06);

    public static final JitterVariableParams ULTIMA_WGS_DN = new JitterVariableParams(
            0.19, 0.195, 0.20, 0.01);

    public static final JitterVariableParams ULTIMA_WGS_3P = new JitterVariableParams(
            0.28, 0.29, 0.30, 0.01);

    // SBX defaults
    public static final JitterVariableParams SBX_WGS_HP = new JitterVariableParams(
            0.14, 0.155, 0.17, 0.02);

    public static final JitterVariableParams SBX_WGS_DN = new JitterVariableParams(
            0.11, 0.15, 0.19, 0.02);

    public static final JitterVariableParams SBX_WGS_3P = new JitterVariableParams(
            0.14, 0.18, 0.22, 0.06);

    private static final List<String> HP_REPEAT_UNITS = List.of("A/T", "C/G");
    private static final List<String> DN_REPEAT_UNITS = List.of("AT/TA", "AC/CA/GT/TG", "AG/GA/CT/TC", "CG/GC");
    private static final List<String> THREE_PLUS_REPEAT_UNIT = List.of(REPEAT_UNIT_3_PLUS_LABEL);

    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS_ILLUMINA = Lists.newArrayList();
    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS_ILLUMINA_HIGH_DEPTH = Lists.newArrayList();
    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS_ULTIMA = Lists.newArrayList();
    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS_SBX = Lists.newArrayList();

    private static void addDefaultJitterParams(
            final List<JitterModelParams> collection, final List<String> repeatUnits, final JitterVariableParams variableParams)
    {
        double intercept = variableParams.Optimal6 - 6 * variableParams.Gradient;

        for(String hpRepeat : repeatUnits)
        {
            collection.add(
                    new JitterModelParams(
                            hpRepeat, ConsensusType.NONE, variableParams.Optimal4, variableParams.Optimal5, variableParams.Optimal6,
                            variableParams.Gradient, intercept, 1));
        }
    }

    static
    {
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA, HP_REPEAT_UNITS, ILLUMINA_WGS_HP);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA, DN_REPEAT_UNITS, ILLUMINA_WGS_DN);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA, THREE_PLUS_REPEAT_UNIT, ILLUMINA_WGS_3P);

        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA_HIGH_DEPTH, HP_REPEAT_UNITS, ILLUMINA_HIGH_DEPTH_HP);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA_HIGH_DEPTH, DN_REPEAT_UNITS, ILLUMINA_HIGH_DEPTH_DN);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ILLUMINA_HIGH_DEPTH, THREE_PLUS_REPEAT_UNIT, ILLUMINA_HIGH_DEPTH_3P);

        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ULTIMA, HP_REPEAT_UNITS, ULTIMA_WGS_HP);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ULTIMA, DN_REPEAT_UNITS, ULTIMA_WGS_DN);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_ULTIMA, THREE_PLUS_REPEAT_UNIT, ULTIMA_WGS_3P);

        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_SBX, HP_REPEAT_UNITS, SBX_WGS_HP);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_SBX, DN_REPEAT_UNITS, SBX_WGS_DN);
        addDefaultJitterParams(DEFAULT_JITTER_PARAMS_SBX, THREE_PLUS_REPEAT_UNIT, SBX_WGS_3P);
    }

    // DEPRECATED: to be removed
    public static double DEFAULT_JITTER_SCALE_GRADIENT = 0.04;
    public static double DEFAULT_JITTER_SCALE_4_HP = 0.05;
    public static double DEFAULT_JITTER_SCALE_5_HP = 0.09;
    public static double DEFAULT_JITTER_SCALE_6_HP = 0.13;
    public static double DEFAULT_JITTER_SCALE_INTERCEPT_HP = DEFAULT_JITTER_SCALE_6_HP - 6 * DEFAULT_JITTER_SCALE_GRADIENT;
    public static double DEFAULT_JITTER_SCALE_4_DN = 0.11;
    public static double DEFAULT_JITTER_SCALE_5_DN = 0.15;
    public static double DEFAULT_JITTER_SCALE_6_DN = 0.19;
    public static double DEFAULT_JITTER_SCALE_INTERCEPT_DN = DEFAULT_JITTER_SCALE_6_DN - 6 * DEFAULT_JITTER_SCALE_GRADIENT;
    public static double DEFAULT_JITTER_SCALE_4_3P = 0.15;
    public static double DEFAULT_JITTER_SCALE_5_3P = 0.19;
    public static double DEFAULT_JITTER_SCALE_6_3P = 0.23;
    public static double DEFAULT_JITTER_SCALE_INTERCEPT_3P = DEFAULT_JITTER_SCALE_6_3P - 6 * DEFAULT_JITTER_SCALE_GRADIENT;

    // high depth equivalents
    // public static double DEFAULT_HD_JITTER_SCALE_GRADIENT = 0.06;
    public static double DEFAULT_HD_JITTER_SCALE_4_HP = 0.1;
    public static double DEFAULT_HD_JITTER_SCALE_5_HP = 0.15;
    public static double DEFAULT_HD_JITTER_SCALE_6_HP = 0.2;
    public static double DEFAULT_HD_JITTER_SCALE_INTERCEPT_HP = DEFAULT_HD_JITTER_SCALE_6_HP - 6 * DEFAULT_HD_JITTER_SCALE_GRADIENT;
    public static double DEFAULT_HD_JITTER_SCALE_4_DN = 0.15;
    public static double DEFAULT_HD_JITTER_SCALE_5_DN = 0.2;
    public static double DEFAULT_HD_JITTER_SCALE_6_DN = 0.25;
    public static double DEFAULT_HD_JITTER_SCALE_INTERCEPT_DN = DEFAULT_HD_JITTER_SCALE_6_DN - 6 * DEFAULT_HD_JITTER_SCALE_GRADIENT;
    public static double DEFAULT_HD_JITTER_SCALE_4_3P = 0.2;
    public static double DEFAULT_HD_JITTER_SCALE_5_3P = 0.25;
    public static double DEFAULT_HD_JITTER_SCALE_6_3P = 0.3;
    public static double DEFAULT_HD_JITTER_SCALE_INTERCEPT_3P = DEFAULT_HD_JITTER_SCALE_6_3P - 6 * DEFAULT_HD_JITTER_SCALE_GRADIENT;

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_HP_1 = new JitterModelParams(
            "A/T", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_HP, DEFAULT_JITTER_SCALE_5_HP, DEFAULT_JITTER_SCALE_6_HP, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_HP, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_HP_2 = new JitterModelParams(
            "C/G", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_HP, DEFAULT_JITTER_SCALE_5_HP, DEFAULT_JITTER_SCALE_6_HP, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_HP, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_DN_1 = new JitterModelParams(
            "AT/TA", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_DN, DEFAULT_JITTER_SCALE_5_DN, DEFAULT_JITTER_SCALE_6_DN, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_DN_2 = new JitterModelParams(
            "AC/CA/GT/TG", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_DN, DEFAULT_JITTER_SCALE_5_DN, DEFAULT_JITTER_SCALE_6_DN, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_DN_3 = new JitterModelParams(
            "AG/GA/CT/TC", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_DN, DEFAULT_JITTER_SCALE_5_DN, DEFAULT_JITTER_SCALE_6_DN, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_DN_4 = new JitterModelParams(
            "CG/GC", ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_DN, DEFAULT_JITTER_SCALE_5_DN, DEFAULT_JITTER_SCALE_6_DN, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_3P = new JitterModelParams(
            REPEAT_UNIT_3_PLUS_LABEL, ConsensusType.NONE, DEFAULT_JITTER_SCALE_4_3P, DEFAULT_JITTER_SCALE_5_3P, DEFAULT_JITTER_SCALE_6_3P, DEFAULT_JITTER_SCALE_GRADIENT,
            DEFAULT_JITTER_SCALE_INTERCEPT_3P, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_HP_1 = new JitterModelParams(
            "A/T", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_HP, DEFAULT_HD_JITTER_SCALE_5_HP, DEFAULT_HD_JITTER_SCALE_6_HP, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_HP, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_HP_2 = new JitterModelParams(
            "C/G", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_HP, DEFAULT_HD_JITTER_SCALE_5_HP, DEFAULT_HD_JITTER_SCALE_6_HP, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_HP, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_DN_1 = new JitterModelParams(
            "AT/TA", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_DN, DEFAULT_HD_JITTER_SCALE_5_DN, DEFAULT_HD_JITTER_SCALE_6_DN, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_DN_2 = new JitterModelParams(
            "AC/CA/GT/TG", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_DN, DEFAULT_HD_JITTER_SCALE_5_DN, DEFAULT_HD_JITTER_SCALE_6_DN, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_DN_3 = new JitterModelParams(
            "AG/GA/CT/TC", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_DN, DEFAULT_HD_JITTER_SCALE_5_DN, DEFAULT_HD_JITTER_SCALE_6_DN, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_DN_4 = new JitterModelParams(
            "CG/GC", ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_DN, DEFAULT_HD_JITTER_SCALE_5_DN, DEFAULT_HD_JITTER_SCALE_6_DN, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_DN, 1);

    public static final JitterModelParams DEFAULT_HD_JITTER_PARAMS_3P = new JitterModelParams(
            REPEAT_UNIT_3_PLUS_LABEL, ConsensusType.NONE, DEFAULT_HD_JITTER_SCALE_4_3P, DEFAULT_HD_JITTER_SCALE_5_3P, DEFAULT_HD_JITTER_SCALE_6_3P, DEFAULT_HD_JITTER_SCALE_GRADIENT,
            DEFAULT_HD_JITTER_SCALE_INTERCEPT_3P, 1);

    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS = List.of(
            DEFAULT_JITTER_PARAMS_HP_1, DEFAULT_JITTER_PARAMS_HP_2, DEFAULT_JITTER_PARAMS_DN_1,
            DEFAULT_JITTER_PARAMS_DN_2, DEFAULT_JITTER_PARAMS_DN_3, DEFAULT_JITTER_PARAMS_DN_4, DEFAULT_JITTER_PARAMS_3P);

    public static final List<JitterModelParams> DEFAULT_HD_JITTER_PARAMS = List.of(
            DEFAULT_HD_JITTER_PARAMS_HP_1, DEFAULT_HD_JITTER_PARAMS_HP_2, DEFAULT_HD_JITTER_PARAMS_DN_1,
            DEFAULT_HD_JITTER_PARAMS_DN_2, DEFAULT_HD_JITTER_PARAMS_DN_3, DEFAULT_HD_JITTER_PARAMS_DN_4, DEFAULT_HD_JITTER_PARAMS_3P);

}
