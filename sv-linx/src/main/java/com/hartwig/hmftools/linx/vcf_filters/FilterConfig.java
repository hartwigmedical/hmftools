package com.hartwig.hmftools.linx.vcf_filters;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FilterConfig
{
    public final boolean RequirePass;
    public final boolean LogFiltered;
    public final int QualScoreThreshold;

    private static final String REQUIRE_PASS = "require_pass";
    private static final String LOG_FILTERED = "log_filtered";
    private static final String QUAL_SCORE_THRESHOLD = "qs_threshold";

    public FilterConfig(final CommandLine cmd)
    {
        RequirePass = cmd.hasOption(REQUIRE_PASS);
        LogFiltered = cmd.hasOption(LOG_FILTERED);
        QualScoreThreshold = Integer.parseInt(cmd.getOptionValue(QUAL_SCORE_THRESHOLD, "400"));
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(REQUIRE_PASS, false, "Require variants to have filter = PASS");
        options.addOption(QUAL_SCORE_THRESHOLD, true, "Qual score threshold");
        options.addOption(LOG_FILTERED, false, "Log filtered variants");
    }
}
