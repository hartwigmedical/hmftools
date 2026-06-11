package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.CONTIG_SIDECAR;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.CONTIG_SIDECAR_DESC;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.ENSEMBL_DATA_DIR_DESC;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.EXTEND_SOFTCLIP_TAILS;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.EXTEND_SOFTCLIP_TAILS_DESC;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.RESCUE_VIA_SUPP;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.RESCUE_VIA_SUPP_DESC;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.UNMAP_ABOVE_NH;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.UNMAP_ABOVE_NH_DESC;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.UNMAP_BELOW_MAPQ;
import static com.hartwig.hmftools.redux.splice.SpliceLiftBackConfig.UNMAP_BELOW_MAPQ_DESC;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

// Splice liftback options folded into REDUX. Gated by the master SPLICE_LIFTBACK flag: when absent the
// stage is a no-op and REDUX runs byte-for-byte unchanged. Ref genome, bamtool, threads, output and
// logging are NOT re-registered here -- ReduxConfig.registerConfig already owns them and the stage reads
// them off ReduxConfig. Only the splice-specific options live here. Mirrors the SpliceLiftBackConfig
// fields that aren't already on ReduxConfig (see [[project_redux_migration_scoping]]).
public class SpliceStageConfig
{
    public static final String SPLICE_LIFTBACK = "splice_liftback";
    public static final String SPLICE_LIFTBACK_DESC =
            "Run the transcript-contig liftback pre-pass before dedup (produces a genomic-coord BAM)";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated always-unmap regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "merged into the unmapping regions as high-depth so contained reads are unmapped unconditionally";

    public final boolean Enabled;
    public final String ContigSidecarFile;
    public final String EnsemblDataDir;
    public final boolean RescueViaSupp;
    public final boolean ExtendSoftclipTails;
    public final int UnmapAboveNh;
    public final int UnmapBelowMapq;
    public final String RnaUnmapRegionsFile;

    public SpliceStageConfig(final ConfigBuilder configBuilder)
    {
        Enabled = configBuilder.hasFlag(SPLICE_LIFTBACK);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        EnsemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        RescueViaSupp = configBuilder.hasFlag(RESCUE_VIA_SUPP);
        ExtendSoftclipTails = configBuilder.hasFlag(EXTEND_SOFTCLIP_TAILS);
        UnmapAboveNh = configBuilder.getInteger(UNMAP_ABOVE_NH);
        UnmapBelowMapq = configBuilder.getInteger(UNMAP_BELOW_MAPQ);
        RnaUnmapRegionsFile = configBuilder.getValue(RNA_UNMAP_REGIONS);

        // contig_sidecar / ensembl_data_dir are registered as optional (most REDUX runs don't set the
        // flag); enforce them only when the stage is actually switched on.
        if(Enabled)
        {
            if(ContigSidecarFile == null)
                throw new IllegalArgumentException(SPLICE_LIFTBACK + " requires -" + CONTIG_SIDECAR);
            if(EnsemblDataDir == null)
                throw new IllegalArgumentException(SPLICE_LIFTBACK + " requires -" + ENSEMBL_DATA_DIR);
        }
    }

    public boolean enabled() { return Enabled; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(SPLICE_LIFTBACK, SPLICE_LIFTBACK_DESC);
        configBuilder.addPath(CONTIG_SIDECAR, false, CONTIG_SIDECAR_DESC);
        configBuilder.addPath(ENSEMBL_DATA_DIR, false, ENSEMBL_DATA_DIR_DESC);
        configBuilder.addFlag(RESCUE_VIA_SUPP, RESCUE_VIA_SUPP_DESC);
        configBuilder.addFlag(EXTEND_SOFTCLIP_TAILS, EXTEND_SOFTCLIP_TAILS_DESC);
        configBuilder.addInteger(UNMAP_ABOVE_NH, UNMAP_ABOVE_NH_DESC, 0);
        configBuilder.addInteger(UNMAP_BELOW_MAPQ, UNMAP_BELOW_MAPQ_DESC, 0);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
    }
}
