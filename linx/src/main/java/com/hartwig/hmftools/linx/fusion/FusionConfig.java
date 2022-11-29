package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.linx.LinxConfig.CHECK_FUSIONS;
import static com.hartwig.hmftools.linx.LinxConfig.configPathValid;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionMapper.RNA_FILE_SOURCE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionMapper.RNA_FUSIONS_FILE;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionConfig
{
    public final boolean LogReportableOnly;
    public final boolean LogAllPotentials;
    public final boolean WriteAllVisFusions;

    // dynmamic parameters
    public boolean LogInvalidReasons;
    public boolean AllowExonSkipping;
    public boolean RequirePhaseMatch;
    public boolean RequireUpstreamBiotypes;

    private static final String PRE_GENE_BREAKEND_DISTANCE = "fusion_gene_distance";
    private static final String LOG_REPORTABLE_ONLY = "log_reportable_fusions";
    private static final String LOG_ALL_POTENTIALS = "log_all_potential_fusions";
    public static final String LOG_INVALID_REASONS = "log_invalid_fusions";
    private static final String SKIP_UNPHASED_FUSIONS = "skip_unphased_fusions";
    public static final String WRITE_NEO_EPITOPES = "write_neo_epitopes";
    private static final String WRITE_ALL_VIS_FUSIONS = "write_all_vis_fusions";

    public static boolean validConfig(final CommandLine cmd)
    {
        return configPathValid(cmd, RNA_FUSIONS_FILE) && configPathValid(cmd, KNOWN_FUSIONS_FILE);
    }

    public FusionConfig(final CommandLine cmd)
    {
        if(cmd.hasOption(PRE_GENE_BREAKEND_DISTANCE))
        {
            PRE_GENE_PROMOTOR_DISTANCE = Integer.parseInt(cmd.getOptionValue(PRE_GENE_BREAKEND_DISTANCE));
        }

        LogReportableOnly = cmd.hasOption(LOG_REPORTABLE_ONLY);
        RequirePhaseMatch = cmd.hasOption(SKIP_UNPHASED_FUSIONS);
        LogAllPotentials = cmd.hasOption(LOG_ALL_POTENTIALS);
        WriteAllVisFusions = cmd.hasOption(WRITE_ALL_VIS_FUSIONS);
        RequireUpstreamBiotypes = true;
        LogInvalidReasons = cmd.hasOption(LOG_INVALID_REASONS);
        AllowExonSkipping = true;
    }

    public FusionConfig()
    {
        LogAllPotentials = false;
        LogReportableOnly = false;
        WriteAllVisFusions = false;
        LogInvalidReasons = false;
        AllowExonSkipping = true;
        RequirePhaseMatch = false;
        RequireUpstreamBiotypes = true;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(PRE_GENE_BREAKEND_DISTANCE, true, "Distance after to a breakend to consider in a gene");
        options.addOption(SKIP_UNPHASED_FUSIONS, false, "Skip unphased fusions");
        options.addOption(WRITE_NEO_EPITOPES, false, "Search for neo-epitopes from fusions");

        options.addOption(RNA_FUSIONS_FILE, true, "Sample RNA fusion data to match vs Linx fusions");
        options.addOption(RNA_FILE_SOURCE, true, "RNA fusion source: ISOFOX, ARRIBA or STARFUSION");

        options.addOption(LOG_REPORTABLE_ONLY, false, "Only write out reportable fusions");
        options.addOption(LOG_ALL_POTENTIALS, false, "Log all potential fusions");
        options.addOption(LOG_INVALID_REASONS, false, "Log reasons for not making a fusion between transcripts");
        options.addOption(WRITE_ALL_VIS_FUSIONS, false, "Write all fusions including non-reportable for visualiser");
    }
}
