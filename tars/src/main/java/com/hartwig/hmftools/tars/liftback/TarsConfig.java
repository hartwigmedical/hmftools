package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
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
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryConfig;

import org.jetbrains.annotations.Nullable;

public record TarsConfig(
        String sampleId,
        List<String> inputBams,
        String refGenomeFile,
        String contigSidecarFile,
        SupplementaryConfig supplementary,
        @Nullable String rnaUnmapRegionsFile,
        String outputDir,
        @Nullable String outputId,
        String bamToolPath,
        int threads
)
{
    // output file id, e.g. <sample>.tars.bam (mirrors Redux's FILE_ID)
    private static final String FILE_ID = "tars";

    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC =
            "Input BAM aligned against ref + transcript-contig FASTA, separated by ',' if multiple";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder";

    // Supplementary merge always runs; these expose its implied-intron bounds for the FP/yield tradeoff.
    public static final String SUPP_IMPLIED_MIN_INTRON_LENGTH = "supp_implied_min_intron_length";
    public static final String SUPP_IMPLIED_MAX_INTRON_LENGTH = "supp_implied_max_intron_length";

    public static final String RNA_UNMAP_REGIONS = "rna_unmap_regions";
    public static final String RNA_UNMAP_REGIONS_DESC =
            "Curated excluded regions (Chromosome/PosStart/PosEnd), e.g. RNA rRNA / 7SL / multi-map zones; "
                    + "a read lifting into any region is excluded post-lift (primary unmapped, supplementaries dropped)";

    public static TarsConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String outputDir = parseOutputDir(configBuilder);
        String bamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        if(outputDir == null || bamToolPath == null)
        {
            TARS_LOGGER.error("missing config: outputDir({}) bamTool({})", outputDir != null, bamToolPath != null);
            System.exit(1);
        }

        if(!checkCreateOutputDir(outputDir))
        {
            TARS_LOGGER.error("failed to create output directory: {}", outputDir);
            System.exit(1);
        }

        return new TarsConfig(
                configBuilder.getValue(SAMPLE),
                Arrays.stream(configBuilder.getValue(INPUT_BAM).split(CONFIG_FILE_DELIM, -1)).collect(Collectors.toList()),
                configBuilder.getValue(REF_GENOME),
                configBuilder.getValue(CONTIG_SIDECAR),
                new SupplementaryConfig(
                        configBuilder.getInteger(SUPP_IMPLIED_MIN_INTRON_LENGTH),
                        configBuilder.getInteger(SUPP_IMPLIED_MAX_INTRON_LENGTH),
                        MAX_SUPP_MERGES,
                        false,
                        MAX_SUPP_READ_OVERLAP),
                configBuilder.getValue(RNA_UNMAP_REGIONS),
                outputDir,
                configBuilder.getValue(OUTPUT_ID),
                bamToolPath,
                parseThreads(configBuilder));
    }

    // common output stem; also namespaces the per-worker shard files so repeated runs into one dir don't clobber each other
    public String filePrefix()
    {
        String prefix = outputDir + sampleId + "." + FILE_ID;
        if(outputId != null)
        {
            prefix += "." + outputId;
        }
        return prefix;
    }

    // <sample>.tars[.<outputId>][.<stage>].bam
    public String formBamFilename(final String stage)
    {
        return filePrefix() + (stage != null ? "." + stage : "") + BAM_EXTENSION;
    }

    public String formUnsortedBam() { return formBamFilename("unsorted"); }

    public String formOutputBam() { return formBamFilename(null); }

    public String formSummaryFile() { return filePrefix() + ".summary" + TSV_EXTENSION; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPaths(INPUT_BAM, true, INPUT_BAM_DESC);
        addRefGenomeFile(configBuilder, true);
        // liftback works off the sidecar, not the version, but existing command lines pass it
        addRefGenomeVersion(configBuilder);
        configBuilder.addPath(CONTIG_SIDECAR, true, CONTIG_SIDECAR_DESC);
        configBuilder.addInteger(
                SUPP_IMPLIED_MIN_INTRON_LENGTH,
                "Supplementary merge: min implied intron length", MIN_IMPLIED_INTRON_LENGTH);
        configBuilder.addInteger(
                SUPP_IMPLIED_MAX_INTRON_LENGTH,
                "Supplementary merge: max implied intron length", MAX_IMPLIED_INTRON_LENGTH);
        configBuilder.addPath(RNA_UNMAP_REGIONS, false, RNA_UNMAP_REGIONS_DESC);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
