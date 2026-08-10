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
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_MERGES;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_READ_OVERLAP;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_IMPLIED_INTRON_LENGTH;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryConfig;

public class TarsConfig
{
    // output file id, e.g. <sample>.tars.bam (mirrors Redux's FILE_ID)
    private static final String FILE_ID = "tars";

    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned against ref + transcript-contig FASTA";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder";

    // Supplementary merge always runs; these expose its implied-intron bounds for the FP/yield tradeoff.
    public static final String SUPP_IMPLIED_MIN_INTRON_LENGTH = "supp_implied_min_intron_length";
    public static final String SUPP_IMPLIED_MAX_INTRON_LENGTH = "supp_implied_max_intron_length";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated excluded regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "a read lifting into any region is excluded post-lift (primary unmapped, supplementaries dropped)";

    public final String SampleId;
    public final String InputBam;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String ContigSidecarFile;
    public final SupplementaryConfig Supplementary;
    public final String RnaUnmapRegionsFile;
    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;
    public final int Threads;

    public TarsConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        Supplementary = new SupplementaryConfig(
                configBuilder.getInteger(SUPP_IMPLIED_MIN_INTRON_LENGTH),
                configBuilder.getInteger(SUPP_IMPLIED_MAX_INTRON_LENGTH),
                MAX_SUPP_MERGES,
                false,
                MAX_SUPP_READ_OVERLAP);
        RnaUnmapRegionsFile = configBuilder.getValue(RNA_UNMAP_REGIONS);
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

    // Common output stem; also namespaces per-worker shard intermediates so repeated runs into one dir don't clobber each other.
    public String filePrefix()
    {
        String prefix = OutputDir + SampleId + "." + FILE_ID;
        if(OutputId != null)
        {
            prefix += "." + OutputId;
        }
        return prefix;
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

    public String formSummaryFile() { return filePrefix() + ".summary" + TSV_EXTENSION; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(INPUT_BAM, true, INPUT_BAM_DESC);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, true, CONTIG_SIDECAR_DESC);
        configBuilder.addInteger(
                SUPP_IMPLIED_MIN_INTRON_LENGTH,
                "Supplementary resolve: min implied intron length", MIN_IMPLIED_INTRON_LENGTH);
        configBuilder.addInteger(
                SUPP_IMPLIED_MAX_INTRON_LENGTH,
                "Supplementary resolve: max implied intron length", MAX_IMPLIED_INTRON_LENGTH);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
