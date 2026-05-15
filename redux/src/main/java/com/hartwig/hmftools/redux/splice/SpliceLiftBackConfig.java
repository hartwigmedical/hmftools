package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM (bwa-mem2 with transcript contigs, or generic genomic BAM e.g. STAR)";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC =
            "Contig sidecar TSV produced by SpliceFastaBuilder. Omit for pass-through mode on a generic genomic BAM.";

    public static final String FILTER_RRNA = "filter_rrna";
    public static final String FILTER_RRNA_DESC =
            "Drop read pairs whose post-lift coords overlap an rRNA/Mt_rRNA gene region. Requires -" + ENSEMBL_DATA_DIR;

    public static final String DEFAULT_OUTPUT_PREFIX = "splice_lifted";
    public static final String TSV_A_SUFFIX = ".liftback.records.tsv";
    public static final String TSV_B_SUFFIX = ".liftback.alignments.tsv";
    public static final String SUMMARY_SUFFIX = ".liftback.summary.tsv";

    public final String InputBam;
    public final String RefGenomeFile;
    public final String ContigSidecarFile;
    public final String EnsemblDataDir;
    public final boolean FilterRrna;
    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;
    public final int Threads;

    public SpliceLiftBackConfig(final ConfigBuilder configBuilder)
    {
        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        EnsemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        FilterRrna = configBuilder.hasFlag(FILTER_RRNA);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        Threads = parseThreads(configBuilder);

        if(OutputDir == null)
            throw new IllegalArgumentException("missing required config: output_dir");

        if(FilterRrna && EnsemblDataDir == null)
            throw new IllegalArgumentException("-" + FILTER_RRNA + " requires -" + ENSEMBL_DATA_DIR);

        if(!checkCreateOutputDir(OutputDir))
            throw new IllegalStateException("failed to create output directory: " + OutputDir);
    }

    public String formUnsortedBam()
    {
        return OutputDir + prefix() + ".unsorted" + BAM_EXTENSION;
    }

    public boolean hasContigSidecar()
    {
        return ContigSidecarFile != null;
    }

    private String prefix()
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
        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, false, CONTIG_SIDECAR_DESC);
        addEnsemblDir(configBuilder, false);
        addRefGenomeVersion(configBuilder);
        configBuilder.addFlag(FILTER_RRNA, FILTER_RRNA_DESC);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
