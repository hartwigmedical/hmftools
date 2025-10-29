package com.hartwig.hmftools.pavereverse;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReversePaveConfig
{
    public final String VcfFile;
    public final String TsvInputFile;
    public final String ServeJsonInputFile;
    public final String TsvOuputFile;
    public final String Mode;
    public final String EnsembleDataDir;
    public RefGenomeSource RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final EnsemblDataCache EnsemblCache;

    public static final String VCF_INPUT_FILE = "vcf_input";
    public static final String TSV_INPUT_FILE = "tsv_input";
    public static final String SERVE_JSON_INPUT_FILE = "serve_json_input";
    public static final String TSV_OUTPUT_FILE = "tsv_output";
    private static final String MODE = "mode";
    static String ROUND_TRIP_MODE = "round_trip";
    static String BATCH_MODE = "batch";
    static String SERVE_JSON_MODE = "serve_json";

    public static final Logger RPV_LOGGER = LogManager.getLogger(ReversePaveConfig.class);

    public ReversePaveConfig(final ConfigBuilder configBuilder)
    {
        Mode = configBuilder.getValue(MODE, BATCH_MODE);
        VcfFile = configBuilder.getValue(VCF_INPUT_FILE);
        TsvInputFile = configBuilder.getValue(TSV_INPUT_FILE);
        TsvOuputFile = configBuilder.getValue(TSV_OUTPUT_FILE);
        ServeJsonInputFile = configBuilder.getValue(SERVE_JSON_INPUT_FILE);
        EnsembleDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        try
        {
            String refGenomeFile = configBuilder.getValue(REF_GENOME);
            IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(refGenomeFile));
            RefGenome = new RefGenomeSource(refGenome);
        }
        catch(IOException e)
        {
            RPV_LOGGER.error("failed to load ref genome: {}", e.toString());
        }
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        EnsemblCache = new EnsemblDataCache(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(MODE, false, "Operation mode");
        configBuilder.addPath(VCF_INPUT_FILE, false, "VCF input file");
        configBuilder.addPath(TSV_INPUT_FILE, false, "TSV input file");
        configBuilder.addPath(SERVE_JSON_INPUT_FILE, false, "Serve json input file");
        configBuilder.addConfigItem(TSV_OUTPUT_FILE, false, "TSV output file");

        addRefGenomeConfig(configBuilder, true);
        addEnsemblDir(configBuilder, true);

        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
