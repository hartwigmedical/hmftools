package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.from;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public record SpliceFastaConfig(
        String ensemblDataDir,
        String refGenomeFile,
        RefGenomeVersion refGenVersion,
        String outputDir,
        @Nullable String outputId)
{
    public static SpliceFastaConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String outputDir = parseOutputDir(configBuilder);

        if(outputDir == null)
        {
            TARS_LOGGER.error("missing config: outputDir");
            System.exit(1);
        }

        if(!checkCreateOutputDir(outputDir))
        {
            TARS_LOGGER.error("failed to create output directory: {}", outputDir);
            System.exit(1);
        }

        return new SpliceFastaConfig(
                configBuilder.getValue(ENSEMBL_DATA_DIR),
                configBuilder.getValue(REF_GENOME),
                configBuilder.hasValue(REF_GENOME_VERSION) ? from(configBuilder) : RefGenomeVersion.V38,
                outputDir,
                configBuilder.getValue(OUTPUT_ID));
    }

    public String formFilename(final String fileId)
    {
        String prefix = outputId != null ? outputId : "ref_genome_v" + refGenVersion.identifier() + "_rna_contigs";
        return outputDir + prefix + fileId;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addEnsemblDir(configBuilder, true);
        addRefGenomeFile(configBuilder, true);
        addRefGenomeVersion(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
