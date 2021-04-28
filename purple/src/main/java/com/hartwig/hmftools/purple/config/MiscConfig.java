package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.purple.config.PurpleConstants.DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE;
import static com.hartwig.hmftools.purple.config.PurpleConstants.DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class MiscConfig
{
    // recovery
    public final int RecoveryMinMateQualScore;
    public final int RecoveryMinSglQualScore;

    public static final String CFG_MIN_MATE_QUAL_SCORE = "recovery_mate_min_qual";
    public static final String CFG_MIN_SGL_QUAL_SCORE = "recovery_sgl_min_qual";

    public MiscConfig(final CommandLine cmd)
    {
        RecoveryMinMateQualScore = getConfigValue(cmd, CFG_MIN_MATE_QUAL_SCORE, DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE);
        RecoveryMinSglQualScore = getConfigValue(cmd, CFG_MIN_SGL_QUAL_SCORE, DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE);
    }

    public static void addOptions(@NotNull Options options)
    {
        options.addOption(
                CFG_MIN_MATE_QUAL_SCORE, true,
                "SV recovery non-SGL min qual score (default " + DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE + ")");

        options.addOption(
                CFG_MIN_SGL_QUAL_SCORE, true,
                "SV recovery SGL min qual score (default " + DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE + ")");
    }

}
