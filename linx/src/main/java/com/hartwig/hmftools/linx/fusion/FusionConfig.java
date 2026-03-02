package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionMapper.RNA_FUSIONS_FILE;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FusionConfig
{
    public final boolean LogReportableOnly;
    public final boolean LogAllPotentials;
    public final boolean WriteAllVisFusions;

    public static boolean LOG_INVALID_REASON = false;

    // dynamic parameters
    public boolean AllowExonSkipping;
    public boolean RequirePhaseMatch;
    public boolean RequireUpstreamBiotypes;

    private static final List<FusionSvGeneDebug> FusionDebugSvGenes = Lists.newArrayList();

    private static final String PRE_GENE_BREAKEND_DISTANCE = "fusion_gene_distance";
    private static final String LOG_REPORTABLE_ONLY = "log_reportable_fusions";
    private static final String LOG_ALL_POTENTIALS = "log_all_potential_fusions";
    public static final String LOG_INVALID_REASONS = "log_invalid_fusions";
    private static final String SKIP_UNPHASED_FUSIONS = "skip_unphased_fusions";
    public static final String WRITE_NEO_EPITOPES = "write_neo_epitopes";
    private static final String WRITE_ALL_VIS_FUSIONS = "write_all_vis_fusions";
    private static final String FUSION_DEBUG_SV_GENES = "fusion_debug_sv_genes";

    private class FusionSvGeneDebug
    {
        public final String GeneName;
        public final int SvId;

        public FusionSvGeneDebug(final String geneName, final int svId)
        {
            GeneName = geneName;
            SvId = svId;
        }
    }

    public FusionConfig(final ConfigBuilder configBuilder)
    {
        PRE_GENE_PROMOTOR_DISTANCE = configBuilder.getInteger(PRE_GENE_BREAKEND_DISTANCE);
        LOG_INVALID_REASON = configBuilder.hasFlag(LOG_INVALID_REASONS);

        LogReportableOnly = configBuilder.hasFlag(LOG_REPORTABLE_ONLY);
        RequirePhaseMatch = configBuilder.hasFlag(SKIP_UNPHASED_FUSIONS);
        LogAllPotentials = configBuilder.hasFlag(LOG_ALL_POTENTIALS);
        WriteAllVisFusions = configBuilder.hasFlag(WRITE_ALL_VIS_FUSIONS);
        RequireUpstreamBiotypes = true;
        AllowExonSkipping = true;

        if(configBuilder.hasValue(FUSION_DEBUG_SV_GENES))
        {
            String[] debugSvGenes = configBuilder.getValue(FUSION_DEBUG_SV_GENES).split(ITEM_DELIM, -1);

            for(String debugSvGenePair : debugSvGenes)
            {
                String[] items = debugSvGenePair.split(":", 2);
                if(items.length == 2)
                {
                    FusionDebugSvGenes.add(new FusionSvGeneDebug(items[0], Integer.parseInt(items[1])));
                }
            }
        }
    }

    public FusionConfig()
    {
        LogAllPotentials = false;
        LogReportableOnly = false;
        WriteAllVisFusions = false;
        AllowExonSkipping = true;
        RequirePhaseMatch = false;
        RequireUpstreamBiotypes = true;
    }

    public static boolean processSvGeneDebug(final String gene, final int svId)
    {
        return FusionDebugSvGenes.isEmpty() || FusionDebugSvGenes.stream().anyMatch(x -> x.GeneName.equals(gene) && x.SvId == svId);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addKnownFusionFileOption(configBuilder);

        configBuilder.addInteger(
                PRE_GENE_BREAKEND_DISTANCE, "Distance after to a breakend to consider in a gene", DEFAULT_PRE_GENE_PROMOTOR_DISTANCE);

        configBuilder.addFlag(SKIP_UNPHASED_FUSIONS, "Skip unphased fusions");
        configBuilder.addFlag(WRITE_NEO_EPITOPES, "Search for neo-epitopes from fusions");

        configBuilder.addConfigItem(RNA_FUSIONS_FILE, "Sample RNA fusion data to match vs Linx fusions");
        configBuilder.addConfigItem(FUSION_DEBUG_SV_GENES, "Debug limit to fusions in these genes");

        configBuilder.addFlag(LOG_REPORTABLE_ONLY, "Only write out reportable fusions");
        configBuilder.addFlag(LOG_ALL_POTENTIALS, "Log all potential fusions");
        configBuilder.addFlag(LOG_INVALID_REASONS, "Log reasons for not making a fusion between transcripts");
        configBuilder.addFlag(WRITE_ALL_VIS_FUSIONS, "Write all fusions including non-reportable for visualiser");
    }
}
