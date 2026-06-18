package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionConfig;

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned against ref + transcript-contig FASTA";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder. Required.";

    // junction rescue (primary + supp merge across junctions) and softclip tail-extension both run always; the
    // following expose their thresholds for tuning the FP/yield tradeoff.
    public static final String RESCUE_MIN_ANCHOR_OVERHANG = "rescue_min_anchor_overhang";
    public static final String RESCUE_MIN_INTRON = "rescue_min_intron";
    public static final String RESCUE_MAX_INTRON = "rescue_max_intron";
    public static final String RESCUE_MAX_CHAIN_DEPTH = "rescue_max_chain_depth";
    public static final String RESCUE_SOFTCLIP_TOLERANCE = "rescue_softclip_tolerance";
    public static final String RESCUE_MAX_BOUNDARY_SHIFT = "rescue_max_boundary_shift";
    public static final String RESCUE_MIN_PARTIAL_MATCH_RUN = "rescue_min_partial_match_run";

    public static final String TAIL_MIN_SOFTCLIP = "tail_min_softclip";
    public static final String TAIL_MIN_EXTENSION = "tail_min_extension";
    public static final String TAIL_MAX_EXTENSION = "tail_max_extension";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated excluded regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "a fragment with a primary in any region is dropped before lifting";

    public static final String WRITE_LIFTBACK_TSV = "write_liftback_tsv";
    public static final String WRITE_LIFTBACK_TSV_DESC =
            "Write per-record liftback debug TSVs (records + alignments). Off by default -- whole-sample TSVs are huge";

    public static final String DEFAULT_OUTPUT_PREFIX = "splice_lifted";
    public static final String TSV_A_SUFFIX = ".liftback.records.tsv";
    public static final String TSV_B_SUFFIX = ".liftback.alignments.tsv";
    public static final String SUMMARY_SUFFIX = ".liftback.summary.tsv";

    public final String InputBam;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String ContigSidecarFile;
    public final RescueConfig Rescue;
    public final TailExtensionConfig TailExtension;
    public final String RnaUnmapRegionsFile;
    public final boolean WriteLiftbackTsv;
    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;
    public final int Threads;

    public SpliceLiftBackConfig(final ConfigBuilder configBuilder)
    {
        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        Rescue = new RescueConfig(
                true,
                configBuilder.getInteger(RESCUE_MIN_INTRON),
                configBuilder.getInteger(RESCUE_MAX_INTRON),
                configBuilder.getInteger(RESCUE_MIN_ANCHOR_OVERHANG),
                configBuilder.getInteger(RESCUE_MAX_CHAIN_DEPTH),
                false,
                configBuilder.getInteger(RESCUE_SOFTCLIP_TOLERANCE),
                configBuilder.getInteger(RESCUE_MAX_BOUNDARY_SHIFT),
                configBuilder.getInteger(RESCUE_MIN_PARTIAL_MATCH_RUN));
        TailExtension = new TailExtensionConfig(
                true,
                configBuilder.getInteger(TAIL_MIN_SOFTCLIP),
                configBuilder.getInteger(TAIL_MIN_EXTENSION),
                configBuilder.getInteger(TAIL_MAX_EXTENSION));
        RnaUnmapRegionsFile = configBuilder.getValue(RNA_UNMAP_REGIONS);
        WriteLiftbackTsv = configBuilder.hasFlag(WRITE_LIFTBACK_TSV);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        Threads = parseThreads(configBuilder);

        if(OutputDir == null)
            throw new IllegalArgumentException("missing required config: output_dir");

        if(!checkCreateOutputDir(OutputDir))
            throw new IllegalStateException("failed to create output directory: " + OutputDir);
    }

    public String formUnsortedBam()
    {
        return OutputDir + prefix() + ".unsorted" + BAM_EXTENSION;
    }

    // the run's output prefix, also used to namespace the per-worker shard intermediates so concurrent or
    // repeated runs into one output dir do not clobber each other's shards.
    public String prefix()
    {
        return OutputId != null ? OutputId : DEFAULT_OUTPUT_PREFIX;
    }

    public String formOutputBam()
    {
        return OutputDir + prefix() + BAM_EXTENSION;
    }

    public String formTsvAFile()
    {
        return OutputDir + prefix() + TSV_A_SUFFIX;
    }

    public String formTsvBFile()
    {
        return OutputDir + prefix() + TSV_B_SUFFIX;
    }

    public String formSummaryFile()
    {
        return OutputDir + prefix() + SUMMARY_SUFFIX;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(INPUT_BAM, true, INPUT_BAM_DESC);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, true, CONTIG_SIDECAR_DESC);
        configBuilder.addInteger(RESCUE_MIN_ANCHOR_OVERHANG,
                "Junction rescue: min anchor overhang each side of a merged junction", RescueConfig.DEFAULT_MIN_ANCHOR_OVERHANG);
        configBuilder.addInteger(RESCUE_MIN_INTRON, "Junction rescue: min intron length", RescueConfig.DEFAULT_MIN_INTRON_LENGTH);
        configBuilder.addInteger(RESCUE_MAX_INTRON, "Junction rescue: max intron length", RescueConfig.DEFAULT_MAX_INTRON_LENGTH);
        configBuilder.addInteger(RESCUE_MAX_CHAIN_DEPTH, "Junction rescue: max supp chain depth", RescueConfig.DEFAULT_MAX_CHAIN_DEPTH);
        configBuilder.addInteger(RESCUE_SOFTCLIP_TOLERANCE,
                "Junction rescue: primary/supp overlap tolerance when snapping to an annotated junction",
                RescueConfig.DEFAULT_SOFTCLIP_TOLERANCE);
        configBuilder.addInteger(RESCUE_MAX_BOUNDARY_SHIFT,
                "Junction rescue: max over-extended bases trimmed when probing an annotated boundary",
                RescueConfig.DEFAULT_MAX_BOUNDARY_SHIFT);
        configBuilder.addInteger(RESCUE_MIN_PARTIAL_MATCH_RUN,
                "Junction rescue: min exon-proximal matched run for partial ref-verify rescue",
                RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        configBuilder.addInteger(TAIL_MIN_SOFTCLIP,
                "Tail extension: min terminal softclip length to consider", TailExtensionConfig.DEFAULT_MIN_SOFTCLIP_LENGTH);
        configBuilder.addInteger(TAIL_MIN_EXTENSION,
                "Tail extension: min ref-matching bases to convert to M", TailExtensionConfig.DEFAULT_MIN_EXTENSION);
        configBuilder.addInteger(TAIL_MAX_EXTENSION,
                "Tail extension: max bases walked into a softclip (caps consuming real junctions)",
                TailExtensionConfig.DEFAULT_MAX_EXTENSION);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
        configBuilder.addFlag(WRITE_LIFTBACK_TSV, WRITE_LIFTBACK_TSV_DESC);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
