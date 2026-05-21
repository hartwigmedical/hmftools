package com.hartwig.hmftools.redux.splice;

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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned against ref + transcript-contig FASTA";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV from SpliceFastaBuilder; omit for pass-through mode";

    public static final String EMIT_SECONDARIES = "emit_secondaries";
    public static final String EMIT_SECONDARIES_DESC = "Input BAM uses 0x100 secondaries (bwa-mem2 -a) instead of XA tags";

    public static final String DROP_LOSING_ALTS = "drop_losing_alts";
    public static final String DROP_LOSING_ALTS_DESC = "Drop secondaries the per-pair discriminator marked Dropped";

    public static final String STAR_MAPQ_LADDER = "star_mapq_ladder";
    public static final String STAR_MAPQ_LADDER_DESC = "Override MAPQ on mapped records using STAR's NH ladder";

    public static final String UNMAP_ABOVE_NH = "unmap_above_nh";
    public static final String UNMAP_ABOVE_NH_DESC = "If > 0, mark records with NH > N as unmapped (STAR default 10)";

    public static final String UNMAP_BELOW_MAPQ = "unmap_below_mapq";
    public static final String UNMAP_BELOW_MAPQ_DESC = "If > 0, mark primary records with final MAPQ < N as unmapped";

    public static final String DEFAULT_OUTPUT_PREFIX = "splice_lifted";
    public static final String TSV_A_SUFFIX = ".liftback.records.tsv";
    public static final String TSV_B_SUFFIX = ".liftback.alignments.tsv";
    public static final String SUMMARY_SUFFIX = ".liftback.summary.tsv";

    public final String InputBam;
    public final String RefGenomeFile;
    public final String ContigSidecarFile;
    public final boolean EmitSecondaries;
    public final boolean DropLosingAlts;
    public final boolean StarMapqLadder;
    public final int UnmapAboveNh;
    public final int UnmapBelowMapq;
    public final String OutputDir;
    public final String OutputId;
    public final String BamToolPath;
    public final int Threads;

    public SpliceLiftBackConfig(final ConfigBuilder configBuilder)
    {
        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        EmitSecondaries = configBuilder.hasFlag(EMIT_SECONDARIES);
        DropLosingAlts = configBuilder.hasFlag(DROP_LOSING_ALTS);
        StarMapqLadder = configBuilder.hasFlag(STAR_MAPQ_LADDER);
        UnmapAboveNh = configBuilder.getInteger(UNMAP_ABOVE_NH);
        UnmapBelowMapq = configBuilder.getInteger(UNMAP_BELOW_MAPQ);
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
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, false, CONTIG_SIDECAR_DESC);
        configBuilder.addFlag(EMIT_SECONDARIES, EMIT_SECONDARIES_DESC);
        configBuilder.addFlag(DROP_LOSING_ALTS, DROP_LOSING_ALTS_DESC);
        configBuilder.addFlag(STAR_MAPQ_LADDER, STAR_MAPQ_LADDER_DESC);
        configBuilder.addInteger(UNMAP_ABOVE_NH, UNMAP_ABOVE_NH_DESC, 0);
        configBuilder.addInteger(UNMAP_BELOW_MAPQ, UNMAP_BELOW_MAPQ_DESC, 0);
        BamToolName.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
