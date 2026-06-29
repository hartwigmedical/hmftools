package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.tars.common.TarsConstants.FILE_ID;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.rescue.RescueConfig;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionConfig;

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned against ref + transcript-contig FASTA";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder";

    // junction rescue (primary + supp merge across junctions) and softclip tail-extension both run always; the
    // following expose their thresholds for tuning the FP/yield tradeoff.
    public static final String RESCUE_MIN_ANCHOR_OVERHANG = "rescue_min_anchor_overhang";
    public static final String RESCUE_MIN_INTRON = "rescue_min_intron";
    public static final String RESCUE_MAX_INTRON = "rescue_max_intron";
    public static final String RESCUE_MAX_CHAIN_DEPTH = "rescue_max_chain_depth";
    public static final String RESCUE_SOFTCLIP_TOLERANCE = "rescue_softclip_tolerance";
    public static final String RESCUE_MAX_BOUNDARY_SHIFT = "rescue_max_boundary_shift";
    public static final String RESCUE_MIN_PARTIAL_MATCH_RUN = "rescue_min_partial_match_run";

    public static final String TAIL_MIN_EXTENSION = "tail_min_extension";
    public static final String TAIL_MAX_EXTENSION = "tail_max_extension";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated excluded regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "a fragment with a primary in any region is dropped before lifting";

    public static final String WRITE_LIFTBACK_TSV = "write_liftback_tsv";
    public static final String WRITE_LIFTBACK_TSV_DESC =
            "Write per-record liftback debug TSVs (records + alignments). Off by default -- whole-sample TSVs are huge";

    // output file-type tokens, e.g. <sample>.tars.summary.tsv
    public static final String SUMMARY_FILE_TYPE = "summary";
    public static final String RECORDS_FILE_TYPE = "records";
    public static final String ALIGNMENTS_FILE_TYPE = "alignments";

    public final String SampleId;
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
        SampleId = configBuilder.getValue(SAMPLE);
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
                configBuilder.getInteger(TAIL_MIN_EXTENSION),
                configBuilder.getInteger(TAIL_MAX_EXTENSION));
        RnaUnmapRegionsFile = configBuilder.getValue(RNA_UNMAP_REGIONS);
        WriteLiftbackTsv = configBuilder.hasFlag(WRITE_LIFTBACK_TSV);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        Threads = parseThreads(configBuilder);

        if(OutputDir == null)
        {
            throw new IllegalArgumentException("missing required config: output_dir");
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            throw new IllegalStateException("failed to create output directory: " + OutputDir);
        }
    }

    // common stem for every output: <dir><sample>.tars[.<outputId>], mirroring Redux. Also namespaces the
    // per-worker shard intermediates so concurrent or repeated runs into one output dir do not clobber each other.
    public String filePrefix()
    {
        String prefix = OutputDir + SampleId + "." + FILE_ID;
        if(OutputId != null)
        {
            prefix += "." + OutputId;
        }
        return prefix;
    }

    // <sample>.tars.<fileType>[.<outputId>].tsv
    public String formFilename(final String fileType)
    {
        String filename = OutputDir + SampleId + "." + FILE_ID + "." + fileType;
        if(OutputId != null)
        {
            filename += "." + OutputId;
        }
        return filename + TSV_EXTENSION;
    }

    // <sample>.tars[.<outputId>][.<stage>].bam
    public String formBamFilename(final String stage)
    {
        String filename = OutputDir + SampleId + "." + FILE_ID;
        if(OutputId != null)
        {
            filename += "." + OutputId;
        }
        if(stage != null)
        {
            filename += "." + stage;
        }
        return filename + BAM_EXTENSION;
    }

    public String formUnsortedBam() { return formBamFilename("unsorted"); }

    public String formOutputBam() { return formBamFilename(null); }

    public String formTsvAFile() { return formFilename(RECORDS_FILE_TYPE); }

    public String formTsvBFile() { return formFilename(ALIGNMENTS_FILE_TYPE); }

    public String formSummaryFile() { return formFilename(SUMMARY_FILE_TYPE); }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(INPUT_BAM, true, INPUT_BAM_DESC);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, true, CONTIG_SIDECAR_DESC);
        configBuilder.addInteger(RESCUE_MIN_ANCHOR_OVERHANG,
                "Junction rescue: min anchor overhang each side of a merged junction", TarsConstants.DEFAULT_MIN_ANCHOR_OVERHANG);
        configBuilder.addInteger(RESCUE_MIN_INTRON, "Junction rescue: min intron length", TarsConstants.DEFAULT_MIN_INTRON_LENGTH);
        configBuilder.addInteger(RESCUE_MAX_INTRON, "Junction rescue: max intron length", TarsConstants.DEFAULT_MAX_INTRON_LENGTH);
        configBuilder.addInteger(RESCUE_MAX_CHAIN_DEPTH, "Junction rescue: max supp chain depth", TarsConstants.DEFAULT_MAX_CHAIN_DEPTH);
        configBuilder.addInteger(RESCUE_SOFTCLIP_TOLERANCE,
                "Junction rescue: primary/supp overlap tolerance when snapping to an annotated junction",
                TarsConstants.DEFAULT_SOFTCLIP_TOLERANCE);
        configBuilder.addInteger(RESCUE_MAX_BOUNDARY_SHIFT,
                "Junction rescue: max over-extended bases trimmed when probing an annotated boundary",
                TarsConstants.DEFAULT_MAX_BOUNDARY_SHIFT);
        configBuilder.addInteger(RESCUE_MIN_PARTIAL_MATCH_RUN,
                "Junction rescue: min exon-proximal matched run for partial ref-verify rescue",
                TarsConstants.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        configBuilder.addInteger(TAIL_MIN_EXTENSION,
                "Tail extension: min softclip length considered and min ref-matching bases converted to M",
                TarsConstants.DEFAULT_MIN_EXTENSION);
        configBuilder.addInteger(TAIL_MAX_EXTENSION,
                "Tail extension: max bases walked into a softclip (caps consuming real junctions)",
                TarsConstants.DEFAULT_MAX_EXTENSION);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
        configBuilder.addFlag(WRITE_LIFTBACK_TSV, WRITE_LIFTBACK_TSV_DESC);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
