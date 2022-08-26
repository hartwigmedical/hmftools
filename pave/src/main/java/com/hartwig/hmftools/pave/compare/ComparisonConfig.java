package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.pave.PaveConstants.ITEM_DELIM;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.ANY;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.CODING_EFFECT;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class ComparisonConfig
{
    public final RefGenomeVersion RefGenVersion;

    public final List<String> SampleIds;
    public final String ReferenceVariantsFile;
    public final boolean OnlyDriverGenes;
    public final boolean OnlyCanonical;
    public final List<ImpactDiffType> ImpactDiffTypes;
    public final boolean WriteTransData;
    public final boolean WriteMatches;
    public final String OutputDir;
    public final String OutputId;
    public final int Threads;

    private static final String REF_VARIANTS_FILE = "ref_variants_file";
    private static final String ONLY_DRIVER_GENES = "only_driver_genes";
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String WRITE_TRANS_DATA = "write_trans_data";
    private static final String WRITE_MATCHES = "write_matches";
    private static final String DIFF_TYPES = "diff_types";

    public ComparisonConfig(final CommandLine cmd)
    {
        SampleIds = loadSampleIdsFile(cmd);

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        ReferenceVariantsFile = cmd.getOptionValue(REF_VARIANTS_FILE);

        ImpactDiffTypes = Lists.newArrayList();

        if(cmd.hasOption(DIFF_TYPES))
        {
            Arrays.stream(cmd.getOptionValue(DIFF_TYPES).split(ITEM_DELIM)).forEach(x -> ImpactDiffTypes.add(ImpactDiffType.valueOf(x)));
        }
        else
        {
            ImpactDiffTypes.add(CODING_EFFECT);
        }

        OnlyCanonical = cmd.hasOption(ONLY_CANONCIAL);
        OnlyDriverGenes = cmd.hasOption(ONLY_DRIVER_GENES);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);
        WriteMatches = cmd.hasOption(WRITE_MATCHES);
        Threads = parseThreads(cmd);
    }

    public boolean checkDiffType(ImpactDiffType type) { return ImpactDiffTypes.contains(type) || ImpactDiffTypes.contains(ANY); }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        addSampleIdFile(options);
        options.addOption(REF_VARIANTS_FILE, true, "File with variants to test against");
        options.addOption(ONLY_DRIVER_GENES, false, "Only compare variants in driver genes");
        options.addOption(ONLY_CANONCIAL, false, "Only compare variants by canonical transcripts");
        addThreadOptions(options);

        options.addOption(
                DIFF_TYPES, true, "Types of diffs to write to file: CODING_EFFECT, EFFECTS, HGVS_CODING, HGVS_PROTEIN, ANY");

        options.addOption(WRITE_TRANS_DATA, false, "Write detailed transcript impact data");
        options.addOption(WRITE_MATCHES, false, "Write matches as well");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version: V37(default) or V38");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        addEnsemblDir(options);
        addDatabaseCmdLineArgs(options);
        addLoggingOptions(options);
        addOutputOptions(options);

        return options;
    }
}
