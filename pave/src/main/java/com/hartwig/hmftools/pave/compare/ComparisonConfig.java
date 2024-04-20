package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.pave.PaveConfig.VCF_FILE;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.ANY;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.CODING_EFFECT;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ComparisonConfig
{
    public final RefGenomeInterface RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final List<String> SampleIds;
    public final String ReferenceVariantsFile;
    public final String SampleVCF;

    public final boolean OnlyDriverGenes;
    public final boolean OnlyCanonical;
    public final List<ImpactDiffType> ImpactDiffTypes;
    public final boolean WriteTransData;
    public final boolean WriteMatches;
    public final String OutputDir;
    public final String OutputId;
    public final int Threads;

    public final GenomeLiftoverCache LiftoverCache;

    private static final String REF_VARIANTS_FILE = "ref_variants_file";
    private static final String LIFTOVER = "liftover";
    private static final String ONLY_DRIVER_GENES = "only_driver_genes";
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String WRITE_TRANS_DATA = "write_trans_data";
    private static final String WRITE_MATCHES = "write_matches";
    private static final String DIFF_TYPES = "diff_types";

    public ComparisonConfig(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
            SampleIds = List.of(configBuilder.getValue(SAMPLE));
        else
            SampleIds = loadSampleIdsFile(configBuilder);

        RefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        RefGenVersion = deriveRefGenomeVersion((RefGenomeSource)RefGenome);

        ReferenceVariantsFile = configBuilder.getValue(REF_VARIANTS_FILE);
        SampleVCF = configBuilder.getValue(VCF_FILE);

        ImpactDiffTypes = Lists.newArrayList();

        if(configBuilder.hasValue(DIFF_TYPES))
        {
            Arrays.stream(configBuilder.getValue(DIFF_TYPES).split(ITEM_DELIM)).forEach(x -> ImpactDiffTypes.add(ImpactDiffType.valueOf(x)));
        }
        else
        {
            ImpactDiffTypes.add(ANY);
        }

        OnlyCanonical = configBuilder.hasFlag(ONLY_CANONCIAL);
        OnlyDriverGenes = configBuilder.hasFlag(ONLY_DRIVER_GENES);

        LiftoverCache = configBuilder.hasFlag(LIFTOVER) ? new GenomeLiftoverCache(true) : null;

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        WriteTransData = configBuilder.hasFlag(WRITE_TRANS_DATA);
        WriteMatches = configBuilder.hasFlag(WRITE_MATCHES);
        Threads = parseThreads(configBuilder);
    }

    public boolean checkDiffType(ImpactDiffType type) { return ImpactDiffTypes.contains(type) || ImpactDiffTypes.contains(ANY); }
    public boolean singleSample() { return SampleIds.size() == 1; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addSampleIdFile(configBuilder, false);

        configBuilder.addPath(REF_VARIANTS_FILE, false, "File with variants to test against");
        configBuilder.addPath(VCF_FILE, false, "Path to sample VCF(s)");

        configBuilder.addFlag(ONLY_DRIVER_GENES, "Only compare variants in driver genes");
        configBuilder.addFlag(ONLY_CANONCIAL, "Only compare variants by canonical transcripts");
        configBuilder.addFlag(LIFTOVER, "Lift-over positions from 37 -> 38");
        addThreadOptions(configBuilder);

        configBuilder.addConfigItem(
                DIFF_TYPES, false,
                "Types of diffs to write to file: " + Arrays.stream(ImpactDiffType.values()).map(x -> x.toString())
                        .collect(Collectors.joining(",")),
                ANY.toString());

        configBuilder.addFlag(WRITE_TRANS_DATA, "Write detailed transcript impact data");
        configBuilder.addFlag(WRITE_MATCHES, "Write matches as well");
        addRefGenomeFile(configBuilder, true);
        configBuilder.addConfigItem(DRIVER_GENE_PANEL_OPTION, DRIVER_GENE_PANEL_OPTION_DESC);
        addEnsemblDir(configBuilder);
        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
    }
}
