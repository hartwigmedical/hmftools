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

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SpliceFastaConfig
{
    public final String EnsemblDataDir;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final String OutputDir;
    public final String OutputId;

    public SpliceFastaConfig(final ConfigBuilder configBuilder)
    {
        EnsemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenVersion = configBuilder.hasValue(REF_GENOME_VERSION) ? from(configBuilder) : RefGenomeVersion.V38;
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(OutputDir == null)
        {
            throw new IllegalArgumentException("missing required config: output_dir");
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            throw new IllegalStateException("failed to create output directory: " + OutputDir);
        }
    }

    public String formFilename(final String fileId)
    {
        String prefix = OutputId != null ? OutputId : "ref_genome_v" + RefGenVersion.identifier() + "_rna_contigs";
        return OutputDir + prefix + fileId;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addEnsemblDir(configBuilder, true);
        addRefGenomeFile(configBuilder, true);
        addRefGenomeVersion(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
