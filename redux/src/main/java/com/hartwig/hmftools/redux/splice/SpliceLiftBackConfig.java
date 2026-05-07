package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SpliceLiftBackConfig
{
    public static final String INPUT_BAM = "input_bam";
    public static final String INPUT_BAM_DESC = "Input BAM aligned to ref + transcript contigs";

    public static final String CONTIG_SIDECAR = "contig_sidecar";
    public static final String CONTIG_SIDECAR_DESC = "Contig sidecar TSV produced by SpliceFastaBuilder";

    public final String InputBam;
    public final String RefGenomeFile;
    public final String ContigSidecarFile;
    public final String OutputDir;
    public final String OutputId;

    private boolean mIsValid;

    public SpliceLiftBackConfig(final ConfigBuilder configBuilder)
    {
        // ConfigBuilder has already enforced presence + existence of every required path before we get here
        mIsValid = true;

        InputBam = configBuilder.getValue(INPUT_BAM);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        ContigSidecarFile = configBuilder.getValue(CONTIG_SIDECAR);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(OutputDir == null)
        {
            RD_LOGGER.error("missing required config: output_dir");
            mIsValid = false;
        }
        else if(!checkCreateOutputDir(OutputDir))
        {
            RD_LOGGER.error("failed to create output directory: {}", OutputDir);
            mIsValid = false;
        }
    }

    public boolean isValid() { return mIsValid; }

    public String formOutputBam()
    {
        String prefix = OutputId != null ? OutputId : "splice_lifted";
        return OutputDir + prefix + BAM_EXTENSION;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(INPUT_BAM, true, INPUT_BAM_DESC);
        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(CONTIG_SIDECAR, true, CONTIG_SIDECAR_DESC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
