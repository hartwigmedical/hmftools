package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.NE_FUSION_COHORT_FILE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class NeoEpitopeAnnotator
{
    private final NeoConfig mConfig;
    private final Map<String,List<NeoEpitopeFusion>> mSampleFusionMap;

    private final EnsemblDataCache mGeneTransCache;
    private final DatabaseAccess mDbAccess;
    private final CohortTpmData mCohortTpmData;

    private BufferedWriter mCohortNeoEpitopeWriter;
    private BufferedWriter mCohortPeptideWriter;

    public NeoEpitopeAnnotator(final CommandLine cmd)
    {
        mConfig = new NeoConfig(cmd);
        mCohortNeoEpitopeWriter = null;
        mCohortPeptideWriter = null;

        mSampleFusionMap = Maps.newHashMap();

        mGeneTransCache = new EnsemblDataCache(cmd.getOptionValue(NeoConfig.GENE_TRANSCRIPTS_DIR), RefGenomeVersion.V37);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        mGeneTransCache.load(false);
        mGeneTransCache.createGeneNameIdMap();

        mDbAccess = DatabaseAccess.createDatabaseAccess(cmd);
        mCohortTpmData = new CohortTpmData(cmd.getOptionValue(NeoConfig.CANCER_TPM_FILE));
    }

    public void run()
    {
        if(mConfig.Samples.isEmpty())
            return;

        if(mConfig.Samples.size() == 1)
        {
            NE_LOGGER.info("processing sample({})", mConfig.Samples.get(0).Id);
        }
        else
        {
            NE_LOGGER.info("processing {} samples", mConfig.Samples.size());

            if(mConfig.WriteCohortFiles)
            {
                mCohortNeoEpitopeWriter = initialiseNeoepitopeWriter(mConfig.OutputDir, null);
                mCohortPeptideWriter = initialisePeptideWriter(mConfig.OutputDir, null);
            }
        }

        // check required inputs and config
        List<NeoSampleTask> sampleTasks = Lists.newArrayList();

        for(final SampleData sample : mConfig.Samples)
        {
            NeoSampleTask sampleTask = new NeoSampleTask(
                    sample, mConfig, mGeneTransCache, mDbAccess, mCohortTpmData, mCohortNeoEpitopeWriter, mCohortPeptideWriter);

            sampleTasks.add(sampleTask);
        }

        if(mConfig.Threads > 1)
        {
            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            sampleTasks.forEach(x -> x.processSample());
        }

        if(mConfig.WriteCohortFiles)
        {
            closeBufferedWriter(mCohortNeoEpitopeWriter);
            closeBufferedWriter(mCohortPeptideWriter);
        }
    }

    private final List<PointMutationData> loadCohortSomaticVariants()
    {
        final List<PointMutationData> pointMutations = Lists.newArrayList();

        if(mConfig.MutationsFile == null)
            return pointMutations;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(mConfig.MutationsFile).toPath());

            if(fileContents.isEmpty())
                return pointMutations;

            final String header = fileContents.get(0);
            fileContents.remove(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            for(String data : fileContents)
            {
                final String[] items = data.split(DELIMITER);

                pointMutations.add(new PointMutationData(
                        items[fieldsIndexMap.get("Chromosome")],
                        Integer.parseInt(items[fieldsIndexMap.get("Position")]),
                        items[fieldsIndexMap.get("Ref")],
                        items[fieldsIndexMap.get("Alt")],
                        items[fieldsIndexMap.get("GeneName")],
                        CodingEffect.valueOf(items[fieldsIndexMap.get("CodingEffect")]),
                        1, -1));
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load common HLA types file({}): {}", mConfig.MutationsFile, e.toString());
        }

        return pointMutations;
    }

    private final List<NeoEpitopeFusion> getSvFusions(final String sampleId)
    {
        if(mSampleFusionMap.isEmpty())
        {
            if(mConfig.SvFusionsDir == null)
                return Lists.newArrayList();

            loadCohortSvNeoEpitopes();
        }

        final List<NeoEpitopeFusion> fusions = mSampleFusionMap.get(sampleId);
        mSampleFusionMap.remove(sampleId);

        return fusions != null ? fusions : Lists.newArrayList();
    }

    private static void populateTpmMedians(
            final CohortTpmData cohortTpmData, final String sampleCancerType,
            final Set<String> upTransNames, final Set<String> downTransNames, final double[] tpmCancer, final double[] tpmCohort)
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            double cancerTotal = 0;
            double cohortTotal = 0;

            for(String transName : fs == FS_UP ? upTransNames : downTransNames)
            {
                final double[] result = cohortTpmData.getTranscriptTpm(transName, sampleCancerType);
                cancerTotal += result[CohortTpmData.CANCER_VALUE];
                cohortTotal += result[CohortTpmData.COHORT_VALUE];
            }

            tpmCancer[fs] = cancerTotal;
            tpmCohort[fs] = cohortTotal;
        }
    }

    private void loadCohortSvNeoEpitopes()
    {
        final String filename = mConfig.SvFusionsDir + NE_FUSION_COHORT_FILE;

        if(Files.exists(Paths.get(filename)))
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                NE_LOGGER.error("empty Linx neo-epitope file({})", filename);
                return;
            }

            int neCount = 0;
            String currentSampleId = "";
            List<NeoEpitopeFusion> fusions = null;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String sampleId = items[0];

                if(!mConfig.Samples.stream().noneMatch(x -> x.Id.equals(sampleId)))
                    continue;

                if(!currentSampleId.equals(sampleId))
                {
                    fusions = Lists.newArrayList();
                    currentSampleId = sampleId;
                    mSampleFusionMap.put(sampleId, fusions);
                }

                fusions.add(NeoEpitopeFusion.fromString(line, true));
                ++neCount;
            }

            NE_LOGGER.debug("loaded {} Linx neo-epitope candidates from file: {}", neCount, filename);
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read Linx neo-epitope file({})", filename, exception.toString());
        }
    }

    public static BufferedWriter initialiseNeoepitopeWriter(final String outputDir, final String sampleId)
    {
        if(outputDir.isEmpty())
            return null;

        try
        {
            String outputFileName = outputDir;

            if(sampleId == null)
                outputFileName += "IMU_NEO_EPITOPES.csv";
            else
                outputFileName += sampleId + NeoConfig.NEO_EPITOPE_FILE_ID;

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(sampleId == null)
                writer.write("SampleId,");

            writer.write(NeoEpitopeFile.header());
            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error initialising neo-epitope output file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeNeoepitopes(
            final BufferedWriter writer, final SampleData sampleData, boolean isCohort, final CohortTpmData cohortTpmData,
            int neId, final NeoEpitope neData, final Set<String> upTransNames, final Set<String> downTransNames)
    {
        if(writer == null)
            return;

        try
        {
            if(isCohort)
                writer.write(String.format("%s,", sampleData.Id));

            final double[] tpmCancer = {0, 0};
            final double[] tpmCohort = {0, 0};
            populateTpmMedians(cohortTpmData, sampleData.CancerType, upTransNames, downTransNames, tpmCancer, tpmCohort);

            final NeoEpitopeFile neFile = neData.toFile(neId, upTransNames, downTransNames, tpmCancer, tpmCohort);
            writer.write(NeoEpitopeFile.toString(neFile));
            writer.newLine();
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    public static BufferedWriter initialisePeptideWriter(final String outputDir, final String sampleId)
    {
        try
        {
            String outputFileName = outputDir;

            if(sampleId == null)
                outputFileName += "IMU_HLA_PEPTIDES.csv";
            else
                outputFileName += sampleId + NeoConfig.HLA_PEPTIDE_FILE_ID;

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(sampleId == null)
                writer.write("SampleId,");

            writer.write("NeId,HlaAllele,Peptide,n_flank,c_flank");
            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error initialising HLA peptide output file: {}", e.toString());
        }

        return null;
    }

    public synchronized static void writePeptideHlaData(
            final BufferedWriter writer, final SampleData sampleData, boolean isCohort,
            final NeoConfig config, int neId, final NeoEpitope neData)
    {
        if(writer == null)
            return;

        if(!config.CommonHlaTypes.isEmpty() && sampleData.HlaTypes.isEmpty())
        {
            sampleData.HlaTypes.addAll(config.CommonHlaTypes);
        }

        if(sampleData.HlaTypes.isEmpty() || config.PeptideLengths[SE_START] == 0 || config.PeptideLengths[SE_END] == 0)
            return;

        try
        {
            if(isCohort)
                writer.write(String.format("%s,", sampleData.Id));

            final List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                    neData.UpstreamAcids, neData.NovelAcid, neData.DownstreamAcids, config.PeptideLengths, config.PeptideFlanks);

            Set<String> uniqueHlaTypes = Sets.newHashSet();

            for(String hlaType : sampleData.HlaTypes)
            {
                if(uniqueHlaTypes.contains(hlaType))
                    continue;

                uniqueHlaTypes.add(hlaType);

                final String predictionHlaType = EpitopeUtils.convertHlaTypeForPredictions(hlaType);

                if(predictionHlaType == null)
                {
                    NE_LOGGER.error("sample({} skipping invalid HLA type: {}", sampleData.Id, hlaType);
                    continue;
                }

                for(PeptideData peptideData : peptides)
                {
                    // skip any peptide which is contained within the upstream wildtype AAs
                    if(neData.UpstreamWildTypeAcids.contains(peptideData.Peptide))
                        continue;

                    // for now skip any upstream peptide containing the 21st AA until the binding prediction routine can handle it
                    if(peptideData.Peptide.contains(AA_SELENOCYSTEINE) || peptideData.UpFlank.contains(AA_SELENOCYSTEINE))
                        continue;

                    writer.write(String.format("%d,%s,%s,%s,%s",
                            neId, predictionHlaType, peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank));
                    writer.newLine();
                }
            }
        }
        catch (final IOException e)
        {
            NE_LOGGER.error("error writing HLA peptide output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NeoEpitopeAnnotator neoEpitopeAnnotator = new NeoEpitopeAnnotator(cmd);
        neoEpitopeAnnotator.run();

        NE_LOGGER.info("Neo-epitope annotations complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
