package com.hartwig.hmftools.imuno.neo_cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_FILE_ID;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.loadSampleIdsFile;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class NeoEpitopeCohort
{
    private final String mOutputDir;
    private final String mSampleDataDir;
    private final List<String> mSampleIds;

    private BufferedWriter mCohortWriter;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";

    public NeoEpitopeCohort(final CommandLine cmd)
    {
        mSampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), mSampleIds);

        mSampleDataDir = cmd.getOptionValue(SAMPLE_DATA_DIR);
        mOutputDir = parseOutputDir(cmd);
        mCohortWriter = null;
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
            return;

        initialiseWriter();

        IM_LOGGER.info("processing {} samples", mSampleIds.size());

        // check required inputs and config
        int processed = 0;

        for(final String sampleId : mSampleIds)
        {
            processSampleNeoEpitopes(sampleId);
            ++processed;

            if(processed > 0 && (processed % 100) == 0)
            {
                IM_LOGGER.info("processed {} samples", processed);
            }
        }

        closeBufferedWriter(mCohortWriter);
    }

    private void processSampleNeoEpitopes(final String sampleId)
    {
        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(mSampleDataDir, sampleId);
            String rnaNeoEpitopeFile = neoEpitopeFile.replace(IM_FILE_ID, ISF_FILE_ID);

            boolean hasRnaData = Files.exists(Paths.get(rnaNeoEpitopeFile));

            final List<String> lines = hasRnaData ?
                    Files.readAllLines(new File(rnaNeoEpitopeFile).toPath()) : Files.readAllLines(new File(neoEpitopeFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            Integer fragIndex = fieldsIndexMap.get("FragmentsNovel");
            Integer baseDepthUp = fieldsIndexMap.get("BaseDepthUp");
            Integer baseDepthDown = fieldsIndexMap.get("BaseDepthDown");

            for(String line : lines)
            {
                NeoEpitopeFile neoData = null;

                try
                {
                    neoData = NeoEpitopeFile.fromString(line, false);
                }
                catch(Exception e)
                {
                    IM_LOGGER.error("sample({}) error({}) parsing neo data: {}", sampleId, e.toString(), line);
                    continue;
                }

                int rnaFragCount = 0;
                int[] rnaBaseDepth = {0, 0};

                if(hasRnaData)
                {
                    final String[] items = line.split(DELIMITER, -1);
                    rnaFragCount = Integer.parseInt(items[fragIndex]);
                    rnaBaseDepth[SE_START] =  Integer.parseInt(items[baseDepthUp]);
                    rnaBaseDepth[SE_END] =  Integer.parseInt(items[baseDepthDown]);
                }

                writeData(sampleId, neoData, hasRnaData, rnaFragCount, rnaBaseDepth);
            }

            IM_LOGGER.debug("sample({}) loaded {} neo-epitopes {}", sampleId, lines.size(), hasRnaData ? "with RNA" : "");
        }
        catch(IOException exception)
        {
            IM_LOGGER.error("failed to read sample({}) neo-epitope file: {}", sampleId, exception.toString());
        }
    }


    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "IMU_NEO_EPITOPES.csv";

            mCohortWriter = createBufferedWriter(outputFileName, false);
            mCohortWriter.write("SampleId,NeId,VariantType,VariantInfo,JunctionCopyNumber");
            mCohortWriter.write(",GeneIdUp,GeneIdDown,GeneNameUp,GeneNameDown");
            mCohortWriter.write(",NmdMin,NmdMax,CodingBasesLengthMin,CodingBasesLengthMax,FusedIntronLength");
            mCohortWriter.write(",SkippedDonors,SkippedAcceptors,UpTranscripts,DownTranscripts");
            mCohortWriter.write(",TpmCancerUp,TpmCohortUp,TpmCancerDown,TpmCohortDown");
            mCohortWriter.write(",HasRna,RnaFragCount,RnaDepthUp,RnaDepthDown");
            mCohortWriter.newLine();
        }
        catch (IOException e)
        {
            IM_LOGGER.error("failed to create neo-epitope writer: {}", e.toString());
        }
    }

    private void writeData(final String sampleId, final NeoEpitopeFile neData, boolean hasRna, int rnaFragCount, final int[] rnaBaseDepth)
    {
        try
        {
            mCohortWriter.write(String.format("%s,%d,%s,%s,%.1f",
                    sampleId, neData.Id, neData.VariantType, neData.VariantInfo, neData.CopyNumber));

            mCohortWriter.write(String.format(",%s,%s,%s,%s",
                    neData.GeneIds[FS_UP], neData.GeneIds[FS_DOWN], neData.GeneNames[FS_UP], neData.GeneNames[FS_DOWN]));

            mCohortWriter.write(String.format(",%d,%d,%d,%d,%d",
                    neData.NmdBases[0], neData.NmdBases[1], neData.CodingBasesLength[0], neData.CodingBasesLength[1], neData.FusedIntronLength));

            mCohortWriter.write(String.format(",%d,%d,%s,%s",
                    neData.SkippedAcceptorsDonors[FS_UP], neData.SkippedAcceptorsDonors[FS_DOWN],
                    neData.Transcripts[FS_UP], neData.Transcripts[FS_DOWN]));

            mCohortWriter.write(String.format(",%6.3e,%6.3e,%6.3e,%6.3e",
                    neData.CancerTpmTotal[FS_UP], neData.CohortTpmTotal[FS_UP], neData.CancerTpmTotal[FS_DOWN], neData.CohortTpmTotal[FS_DOWN]));

            mCohortWriter.write(String.format(",%s,%d,%d,%d",
                    hasRna, rnaFragCount, rnaBaseDepth[SE_START], rnaBaseDepth[SE_END]));

            mCohortWriter.newLine();
        }
        catch (IOException e)
        {
            IM_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }

    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(SAMPLE_DATA_DIR, true, "Directory for sample eno-epitope files");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoEpitopeCohort.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        NeoEpitopeCohort neoEpitopeCohort = new NeoEpitopeCohort(cmd);
        neoEpitopeCohort.run();

        IM_LOGGER.info("Neo-epitope annotations complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
