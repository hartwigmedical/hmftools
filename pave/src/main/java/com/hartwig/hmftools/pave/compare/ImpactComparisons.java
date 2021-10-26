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

            sampleTasks.add(sampleTask);

            sampleTask.call();
        }

        mWriter.close();

        int totalComparisons = sampleTasks.stream().mapToInt(x -> x.totalComparisons()).sum();
        int matchedCount = sampleTasks.stream().mapToInt(x -> x.matchedCount()).sum();

        PV_LOGGER.info("samples({}) total comparisons({}) matched({}) diffs({})",
                mConfig.SampleIds.size(), totalComparisons, matchedCount, totalComparisons - matchedCount);

        // mPerfCounter.logStats();

        PV_LOGGER.info("Pave impact comparison complete");
    }

    /*
    private void loadSampleDatabaseRecords(final String sampleId)
    {
        mTotalComparisons = 0;
        mMatchedCount = 0;

        // retrieve genic positions

        Result<Record20<String, String, Integer, String, String, String, String, String, String, String, Integer, String, String, String,
                String, Integer, Integer, Integer, Byte, String>>
                result = mDbAccess.context()
                .select(SOMATICVARIANT.GENE, SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION,
                        SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.TYPE, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.CANONICALEFFECT, SOMATICVARIANT.CANONICALCODINGEFFECT, SOMATICVARIANT.WORSTCODINGEFFECT,
                        SOMATICVARIANT.GENESEFFECTED, SOMATICVARIANT.CANONICALHGVSCODINGIMPACT, SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT,
                        SOMATICVARIANT.MICROHOMOLOGY, SOMATICVARIANT.REPEATSEQUENCE, SOMATICVARIANT.REPEATCOUNT,
                        SOMATICVARIANT.PHASEDINFRAMEINDEL, SOMATICVARIANT.LOCALPHASESET, SOMATICVARIANT.REPORTED, SOMATICVARIANT.HOTSPOT)
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .orderBy(SOMATICVARIANT.CHROMOSOME,SOMATICVARIANT.POSITION)
                .fetch();

        for(Record record : result)
        {
            //final SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);

            if(mConfig.OnlyDriverGenes)
            {
                String gene = record.getValue(SOMATICVARIANT.GENE);

                if(!mGeneDataCache.isDriverPanelGene(gene))
                    continue;
            }

            processVariant(sampleId, RefVariantData.fromRecord(record));
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET);
        mImpactClassifier.phasedVariants().clear();

        PV_LOGGER.debug("sample({}) processed {} variants, matches({})",
                sampleId, mTotalComparisons, mMatchedCount);
    }

    private void processRefVariantFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final String fileDelim = "\t";
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

            int sampleIndex = fieldsIndexMap.get("sampleId");

            int chrIndex = fieldsIndexMap.get("chromosome");
            int posIndex = fieldsIndexMap.get("position");
            int refIndex = fieldsIndexMap.get("ref");
            int altIndex = fieldsIndexMap.get("alt");
            int typeIndex = fieldsIndexMap.get("type");
            int geneIndex = fieldsIndexMap.get("gene");
            int canonicalEffectIndex = fieldsIndexMap.get("canonicalEffect");
            int canonicalCodingEffectIndex = fieldsIndexMap.get("canonicalCodingEffect");
            int worstCodingEffectIndex = fieldsIndexMap.get("worstCodingEffect");
            int genesAffectedIndex = fieldsIndexMap.get("genesEffected");
            int canonicalHgvsCodingImpactIndex = fieldsIndexMap.get("canonicalHgvsCodingImpact");
            int canonicalHgvsProteinImpactIndex = fieldsIndexMap.get("canonicalHgvsProteinImpact");
            int microhomologyIndex = fieldsIndexMap.get("microhomology");
            int repeatSequenceIndex = fieldsIndexMap.get("repeatSequence");
            int repeatCountIndex = fieldsIndexMap.get("repeatCount");
            int phasedInframeIndelIndex = fieldsIndexMap.get("phasedInframeIndel");
            int localPhaseSetIndex = fieldsIndexMap.get("localPhaseSet");
            int reportedIndex = fieldsIndexMap.get("reported");
            Integer hotspotIndex = fieldsIndexMap.get("hotspot");

            String currentSample = "";
            int samplesProcessed = 0;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(fileDelim, -1);

                String sampleId = items[sampleIndex];

                if(!currentSample.equals(sampleId))
                {
                    if(!mConfig.SampleIds.isEmpty() && !mConfig.SampleIds.contains(sampleId))
                    {
                        if(samplesProcessed == mConfig.SampleIds.size())
                        {
                            PV_LOGGER.info("all {} samples processed", samplesProcessed);
                            break;
                        }
                        else
                        {
                            continue;
                        }
                    }

                    currentSample = sampleId;
                    ++samplesProcessed;

                    if(samplesProcessed > 0 && (samplesProcessed % 100) == 0)
                    {
                        PV_LOGGER.info("processed {} samples", samplesProcessed);
                    }
                }

                int localPhaseSet = !items[localPhaseSetIndex].equals("NULL") ? Integer.parseInt(items[localPhaseSetIndex]) : NO_LOCAL_PHASE_SET;

                RefVariantData variant = new RefVariantData(
                        items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                        VariantType.valueOf(items[typeIndex]), items[geneIndex], items[canonicalEffectIndex],
                        items[canonicalCodingEffectIndex].isEmpty() ? NONE : CodingEffect.valueOf(items[canonicalCodingEffectIndex]),
                        CodingEffect.valueOf(items[worstCodingEffectIndex]), Integer.parseInt(items[genesAffectedIndex]),
                        items[canonicalHgvsCodingImpactIndex], items[canonicalHgvsProteinImpactIndex],
                        items[microhomologyIndex], items[repeatSequenceIndex], Integer.parseInt(items[repeatCountIndex]),
                        Boolean.parseBoolean(items[phasedInframeIndelIndex]),
                        localPhaseSet, items[reportedIndex].equals("1"),
                        hotspotIndex != null ? items[hotspotIndex].equals(Hotspot.HOTSPOT.toString()) : false);

                processVariant(sampleId, variant);
            }

            processPhasedVariants(NO_LOCAL_PHASE_SET);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read ref variant data file: {}", e.toString());
        }
    }
    */

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
