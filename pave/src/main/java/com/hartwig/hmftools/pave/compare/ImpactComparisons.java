package com.hartwig.hmftools.pave.compare;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;
import static com.hartwig.hmftools.pave.PaveApplication.findVariantImpacts;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasCodingEffectDiff;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.CODING_EFFECT;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_CODING;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.HGVS_PROTEIN;
import static com.hartwig.hmftools.pave.compare.ImpactDiffType.REPORTED;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantDAO;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.ImpactClassifier;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.VariantImpactBuilder;
import com.hartwig.hmftools.pave.VariantTransImpact;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Record20;
import org.jooq.Result;

public class ImpactComparisons
{
    private final ComparisonConfig mConfig;
    private final GeneDataCache mGeneDataCache;
    private final DatabaseAccess mDbAccess;
    private final ComparisonWriter mWriter;
    private final RefGenomeInterface mRefGenome;

    public ImpactComparisons(final CommandLine cmd)
    {
        mConfig = new ComparisonConfig(cmd);

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion, cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION),
                true, false);

        mRefGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));
        mDbAccess = createDatabaseAccess(cmd);

        mWriter = new ComparisonWriter(mGeneDataCache, mConfig);
    }

    public void run()
    {
        if(mConfig.SampleIds.isEmpty() && mConfig.ReferenceVariantsFile == null)
        {
            PV_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        if(mDbAccess == null && mConfig.ReferenceVariantsFile == null)
        {
            PV_LOGGER.error("neither DB nor ref variants file configured, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache(mConfig.OnlyCanonical, mConfig.OnlyDriverGenes))
        {
            PV_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        Map<String,List<RefVariantData>> sampleVariantsCache = DataLoader.processRefVariantFile(mConfig.ReferenceVariantsFile);

        if(mConfig.SampleIds.isEmpty() && !sampleVariantsCache.isEmpty())
        {
            sampleVariantsCache.keySet().forEach(x -> mConfig.SampleIds.add(x));
        }

        List<SampleComparisonTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new SampleComparisonTask(
                        i, mConfig, mRefGenome, mDbAccess, mWriter, mGeneDataCache, sampleVariantsCache));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            SampleComparisonTask sampleTask = new SampleComparisonTask(
                    0, mConfig, mRefGenome, mDbAccess, mWriter, mGeneDataCache, sampleVariantsCache);

            sampleTask.getSampleIds().addAll(mConfig.SampleIds);

            sampleTasks.add(sampleTask);

            sampleTask.call();
        }

        mWriter.close();

        int totalComparisons = sampleTasks.stream().mapToInt(x -> x.totalComparisons()).sum();
        int matchedCount = sampleTasks.stream().mapToInt(x -> x.matchedCount()).sum();

        PV_LOGGER.info("samples({}) total comparisons({}) matched({}) diffs({})",
                mConfig.SampleIds.size(), totalComparisons, matchedCount, totalComparisons - matchedCount);

        if(PV_LOGGER.isDebugEnabled())
        {
            Map<String,PerformanceCounter> combinedPerfCounters = sampleTasks.get(0).getPerfCounters();

            for(int i = 1; i < sampleTasks.size(); ++i)
            {
                Map<String,PerformanceCounter> saPerfCounters = sampleTasks.get(i).getPerfCounters();

                for(Map.Entry<String,PerformanceCounter> entry : combinedPerfCounters.entrySet())
                {
                    PerformanceCounter combinedPc = entry.getValue();
                    PerformanceCounter saPc = saPerfCounters.get(entry.getKey());

                    if(combinedPc != null)
                        combinedPc.merge(saPc);
                }
            }

            for(Map.Entry<String,PerformanceCounter> entry : combinedPerfCounters.entrySet())
            {
                entry.getValue().logStats();
            }
        }

        PV_LOGGER.info("Pave impact comparison complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = ComparisonConfig.createOptions();

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ImpactComparisons impactComparison = new ImpactComparisons(cmd);
        impactComparison.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
