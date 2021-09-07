package com.hartwig.hmftools.svtools.cohort;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.svtools.cohort.LineInsertSiteData.INSERT_TYPE_PSD;
import static com.hartwig.hmftools.svtools.cohort.LineInsertSiteData.PROGRAM_LINX;
import static com.hartwig.hmftools.svtools.cohort.LineInsertSiteData.PROGRAM_PCAWG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class LineInsertSiteCompare
{
    private static final Logger LOGGER = LogManager.getLogger(CohortLineElements.class);

    private final Map<String,List<LineInsertSiteData>> mSampleInsertSiteData;

    private final String mOutputDir;
    private final String mLinxDataFile;
    private final String mExtDataFile;
    private BufferedWriter mWriter;

    public LineInsertSiteCompare(final CommandLine cmd)
    {
        mSampleInsertSiteData = Maps.newHashMap();

        mOutputDir = parseOutputDir(cmd);
        mLinxDataFile = cmd.getOptionValue(LINX_DATA_FILE);
        mExtDataFile = cmd.getOptionValue(EXT_DATA_FILE);
        mWriter = null;
    }

    public void run()
    {
        loadLinxDataFile(mLinxDataFile);
        loadExternalDataFile(mExtDataFile);

        produceResults();
    }

    private void loadLinxDataFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            int lineSvCount = 0;

            for(final String line : fileContents)
            {
                ++lineSvCount;

                LineInsertSiteData insertData = LineInsertSiteData.fromLinx(fieldsIndexMap, line);

                if(insertData != null)
                {
                    List<LineInsertSiteData> lineDataList = mSampleInsertSiteData.get(insertData.SampleId);

                    if(lineDataList == null)
                    {
                        lineDataList = Lists.newArrayList();
                        mSampleInsertSiteData.put(insertData.SampleId, lineDataList);
                    }

                    lineDataList.add(insertData);
                }
            }

            LNX_LOGGER.info("loaded {} linx line chain items from file: {}", lineSvCount, filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read linx line chain data file({})", filename, exception.toString());
        }
    }

    private void loadExternalDataFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            int lineSvCount = 0;

            for(final String line : fileContents)
            {
                ++lineSvCount;

                LineInsertSiteData insertData = LineInsertSiteData.fromExternal(fieldsIndexMap, line);

                if(insertData != null)
                {
                    if(insertData.InsertType.equals(INSERT_TYPE_PSD)) // pseudo-gene insertions skipped from comparison
                        continue;

                    List<LineInsertSiteData> lineDataList = mSampleInsertSiteData.get(insertData.SampleId);

                    if(lineDataList == null)
                    {
                        lineDataList = Lists.newArrayList();
                        mSampleInsertSiteData.put(insertData.SampleId, lineDataList);
                    }

                    lineDataList.add(insertData);
                }
            }

            LNX_LOGGER.info("loaded {} external line insertion sites from file: {}", lineSvCount, filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read external line insertion data file({})", filename, exception.toString());
        }
    }

    private void produceResults()
    {
        for(Map.Entry<String,List<LineInsertSiteData>> entry : mSampleInsertSiteData.entrySet())
        {
            final List<LineInsertSiteData> insertDataList = entry.getValue();

            final List<LineInsertSiteData> linxDataList = insertDataList.stream().filter(x -> x.Program.equals(PROGRAM_LINX)).collect(Collectors.toList());
            final List<LineInsertSiteData> pcawgDataList = insertDataList.stream().filter(x -> x.Program.equals(PROGRAM_PCAWG)).collect(Collectors.toList());

            for(final LineInsertSiteData linxData : linxDataList)
            {
                final LineInsertSiteData pcawgData = pcawgDataList.stream().filter(x -> x.matches(linxData)).findFirst().orElse(null);

                if(pcawgData != null)
                {
                    writeMatchResults(linxData, pcawgData);
                    pcawgDataList.remove(pcawgData);
                }
                else
                {
                    writeMatchResults(linxData, null);
                }
            }

            for(final LineInsertSiteData pcawgData : pcawgDataList)
            {
                writeMatchResults(null, pcawgData);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void writeMatchResults(final LineInsertSiteData linxData, final LineInsertSiteData pcawgData)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFileName = mOutputDir + "LINE_INSERT_COMPARE.csv";

                mWriter = createBufferedWriter(outputFileName, false);

                mWriter.write("SampleId,MatchType,InsertChr,InsertPosStart,InsertPosEnd,InsertType");
                mWriter.write(",SourceChromosome,SourcePosStart,SourcePosEnd,ClusterId,ChainId,ChainDesc,HasInversion");

                mWriter.newLine();
            }

            final LineInsertSiteData refData = linxData != null ? linxData : pcawgData;

            final String matchType = linxData != null && pcawgData != null ? "MATCH" : (linxData != null ? PROGRAM_LINX : PROGRAM_PCAWG);

            mWriter.write(String.format("%s,%s,%s,%d,%d,%s",
                    refData.SampleId, matchType,
                    refData.InsertSite.Chromosome, refData.InsertSite.start(), refData.InsertSite.end(), refData.InsertType));

            final String sourceStr = refData.SourceSite != null ?
                    String.format("%s,%d,%d", refData.SourceSite.Chromosome, refData.SourceSite.start(), refData.SourceSite.end()) : "-1,-1,-1";

            mWriter.write(String.format(",%s,%d,%d,%s,%s",
                    sourceStr, refData.ClusterId, refData.ChainId, refData.ChainDesc, refData.HasInversion));

            mWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("failed to write to line cluster data: {}", e.toString());
        }
    }


    private static final String LINX_DATA_FILE = "linx_data_file";
    private static final String EXT_DATA_FILE = "ext_data_file";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(LINX_DATA_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(EXT_DATA_FILE, true, "External LINE data sample counts");
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        LineInsertSiteCompare lineInsertSiteCompare = new LineInsertSiteCompare(cmd);
        lineInsertSiteCompare.run();

        LOGGER.info("LINE insert site comparision complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
