package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ComparisonConfig
{
    public final RefGenomeVersion RefGenVersion;

    public final List<String> SampleIds;
    public final String ReferenceVariantsFile;
    public final boolean OnlyDriverGenes;
    public final boolean OnlyCanonical;
    public final boolean WriteTransData;
    public final String OutputDir;
    public final String OutputId;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String REF_VARIANTS_FILE = "ref_variants_file";
    private static final String ONLY_DRIVER_GENES = "only_driver_genes";
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String OUTPUT_ID = "output_id";
    private static final String WRITE_TRANS_DATA = "write_trans_data";

    public static final Logger PV_LOGGER = LogManager.getLogger(ComparisonConfig.class);

    public ComparisonConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_ID_FILE))
            SampleIds.addAll(loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE)));

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        ReferenceVariantsFile = cmd.getOptionValue(REF_VARIANTS_FILE);

        OnlyCanonical = cmd.hasOption(ONLY_CANONCIAL);
        OnlyDriverGenes = cmd.hasOption(ONLY_DRIVER_GENES);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE_ID_FILE, true, "Sample ID file");
        options.addOption(REF_VARIANTS_FILE, true, "File with variants to test against");
        options.addOption(ONLY_DRIVER_GENES, false, "Only compare variants in driver genes");
        options.addOption(ONLY_CANONCIAL, false, "Only compare variants by canonical transcripts");
        options.addOption(WRITE_TRANS_DATA, false, "Write detailed transcript impact data");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version: V37(default) or V38");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        addEnsemblDir(options);
        addDatabaseCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file identifier");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }
}
