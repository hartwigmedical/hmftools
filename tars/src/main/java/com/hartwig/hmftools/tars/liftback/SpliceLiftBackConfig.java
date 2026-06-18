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

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned against ref + transcript-contig FASTA";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder. Required.";

    public static final String RESCUE_VIA_SUPP = "rescue_via_supp";
    public static final String RESCUE_VIA_SUPP_DESC =
            "Merge primary + supplementary across annotated junctions when read coverage complements cleanly";

    public static final String EXTEND_SOFTCLIP_TAILS = "extend_softclip_tails";
    public static final String EXTEND_SOFTCLIP_TAILS_DESC =
            "Convert ref-matching bases at terminal softclips into M (recovers bwa tail-trim losses)";

    public static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";
    public static final String ENSEMBL_DATA_DIR_DESC =
            "Ensembl data cache directory (annotated exons + junctions). Required.";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated excluded regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "a fragment with a primary in any region is dropped before lifting";

    public static final String WRITE_LIFTBACK_TSV = "write_liftback_tsv";
    public static final String WRITE_LIFTBACK_TSV_DESC =
            "Write per-record liftback debug TSVs (records + alignments). Off by default -- whole-sample TSVs are huge";

    public static final String SORT_BAMTOOL_PATH = "sort_bamtool";
    public static final String SORT_BAMTOOL_PATH_DESC =
            "Path to sambamba or samtools used for the final sort only; defaults to -" + BAMTOOL_PATH
                    + " (concat stays samtools)";

    public static final String DEFAULT_OUTPUT_PREFIX = "splice_lifted";
    public static final String TSV_A_SUFFIX = ".liftback.records.tsv";
    public static final String TSV_B_SUFFIX = ".liftback.alignments.tsv";
    public static final String SUMMARY_SUFFIX = ".liftback.summary.tsv";

    public final String InputBam;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String ContigSidecarFile;
    public final boolean RescueViaSupp;
    public final boolean ExtendSoftclipTails;
    public final String EnsemblDataDir;
    public final String RnaUnmapRegionsFile;
    public final boolean WriteLiftbackTsv;
    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;
    public final String SortBamToolPath;
    public final int Threads;

    public SpliceLiftBackConfig(final ConfigBuilder configBuilder)
    {
        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        RescueViaSupp = configBuilder.hasFlag(RESCUE_VIA_SUPP);
        ExtendSoftclipTails = configBuilder.hasFlag(EXTEND_SOFTCLIP_TAILS);
        EnsemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        RnaUnmapRegionsFile = configBuilder.getValue(RNA_UNMAP_REGIONS);
        WriteLiftbackTsv = configBuilder.hasFlag(WRITE_LIFTBACK_TSV);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        SortBamToolPath = configBuilder.hasValue(SORT_BAMTOOL_PATH) ? configBuilder.getValue(SORT_BAMTOOL_PATH) : BamToolPath;
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
        configBuilder.addFlag(RESCUE_VIA_SUPP, RESCUE_VIA_SUPP_DESC);
        configBuilder.addFlag(EXTEND_SOFTCLIP_TAILS, EXTEND_SOFTCLIP_TAILS_DESC);
        configBuilder.addPath(ENSEMBL_DATA_DIR, true, ENSEMBL_DATA_DIR_DESC);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
        configBuilder.addFlag(WRITE_LIFTBACK_TSV, WRITE_LIFTBACK_TSV_DESC);
        BamToolName.addConfig(configBuilder);
        configBuilder.addPath(SORT_BAMTOOL_PATH, false, SORT_BAMTOOL_PATH_DESC);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
