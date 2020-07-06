package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.linx.LinxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_SPLIT;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.LINE_ELEMENT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.fromString;
import static com.hartwig.hmftools.svtools.germline.GermlineVcfConfig.LOG_DEBUG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

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

public class CohortLineElements
{
    private static final Logger LOGGER = LogManager.getLogger(CohortLineElements.class);

    private final Map<String,Map<Integer,LineClusterData>> mSampleClusterLineData;
    private final Map<SvRegion,Integer> mExtLineSampleCounts;

    private final String mOutputDir;
    private final String mSvDataFile;
    private final String mExtDataFile;

    public CohortLineElements(final CommandLine cmd)
    {
        mSampleClusterLineData = Maps.newHashMap();
        mExtLineSampleCounts = Maps.newHashMap();

        mOutputDir = cmd.getOptionValue(DATA_OUTPUT_DIR);
        mSvDataFile = cmd.getOptionValue(SV_DATA_FILE);
        mExtDataFile = cmd.getOptionValue(EXT_DATA_FILE);
    }

    public void run()
    {
        if(mSvDataFile == null)
            return;

        loadLineElementsFile(mSvDataFile);
        loadExternalLineDataFile(mExtDataFile);
        produceResults();
    }

    private void loadLineElementsFile(final String filename)
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

                final String[] items = line.split(",");

                final String sampleId = items[fieldsIndexMap.get("SampleId")];
                final int clusterId = Integer.parseInt(items[fieldsIndexMap.get("ClusterId")]);

                final String[] chromosomes = new String[] { items[fieldsIndexMap.get("ChrStart")], items[fieldsIndexMap.get("ChrEnd")] };

                final int[] positions =
                        new int[] { Integer.parseInt(items[fieldsIndexMap.get("PosStart")]), Integer.parseInt(items[fieldsIndexMap.get("PosEnd")]) };

                final LineElementType[] lineTypes
                        = new LineElementType[] { fromString(items[fieldsIndexMap.get("LEStart")]), fromString(items[fieldsIndexMap.get("LEEnd")]) };

                processLineSv(sampleId, clusterId, chromosomes, positions, lineTypes);
            }

            LNX_LOGGER.info("loaded {} known line elements from file: {}", lineSvCount, filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    private void loadExternalLineDataFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileContents)
            {
                final String[] items = line.split(",");

                final String chromosome = items[fieldsIndexMap.get("Chromosome")];

                final int[] positions =
                        new int[] { Integer.parseInt(items[fieldsIndexMap.get("PosStart")]), Integer.parseInt(items[fieldsIndexMap.get("PosEnd")]) };

                final int sampleCount = Integer.parseInt(items[fieldsIndexMap.get("PcawgSampleCount")]);

                mExtLineSampleCounts.put(new SvRegion(chromosome, positions), sampleCount);
            }

            LNX_LOGGER.info("loaded {} external line data items from file: {}", mExtLineSampleCounts.size(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read line element CSV file({})", filename);
        }

    }

    private void processLineSv(
            final String sampleId, final int clusterId, final String[] chromosomes, final int[] positions, final LineElementType[] lineTypes)
    {
        Map<Integer,LineClusterData> sampleData = mSampleClusterLineData.get(sampleId);

        if(sampleData == null)
        {
            sampleData = Maps.newHashMap();
            mSampleClusterLineData.put(sampleId, sampleData);
        }

        LineClusterData clusterData = sampleData.get(clusterId);

        if(clusterData != null)
        {
            clusterData.addSvData(chromosomes, positions, lineTypes);
            return;
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(chromosomes[se].equals("0"))
                continue;

            if(lineTypes[se] != LineElementType.NONE)
            {
                if(clusterData == null)
                {
                    clusterData = new LineClusterData(
                            sampleId, clusterId, new SvRegion(chromosomes[se], positions[se], positions[se]), lineTypes[se]);
                }
                else
                {
                    clusterData.addSourceData(chromosomes[se], positions[se], lineTypes[se]);
                }
            }
        }

        if(clusterData == null)
        {
            LOGGER.warn("sample({}) cluster({}) has no source elements", sampleId, clusterId);
            return;
        }

        sampleData.put(clusterId, clusterData);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(chromosomes[se].equals("0"))
                continue;

            if(lineTypes[se] == LineElementType.NONE)
            {
                clusterData.addInsertData(chromosomes[se], positions[se]);
            }
        }
    }

    private void produceResults()
    {
        if(mSampleClusterLineData.isEmpty())
            return;

        final List<LineClusterData> combinedLineData = Lists.newArrayList();

        for(Map.Entry<String,Map<Integer,LineClusterData>> entry : mSampleClusterLineData.entrySet())
        {
            final String sampleId = entry.getKey();
            final Map<Integer,LineClusterData> clusterMap = entry.getValue();

            for(LineClusterData clusterData : clusterMap.values())
            {
                LineClusterData matchedRegion = combinedLineData.stream().filter(x -> x.matches(clusterData)).findFirst().orElse(null);

                if(matchedRegion == null)
                {
                    combinedLineData.add(clusterData);
                }
                else
                {
                    matchedRegion.MatchedClusters.add(clusterData);
                }
            }
        }

        writeResults(combinedLineData);
    }

    private void writeResults(final List<LineClusterData> combinedLineData)
    {
        try
        {
            String outputFileName = mOutputDir + "LINE_SOURCE_DATA.csv";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("LineId,Type,Chromosome,PosStart,PosEnd");
            writer.write(",SampleCount,TotalInserts,PcawgSampleCount");
            writer.write(",SampleId,ClusterId,SamplePosStart,SamplePosEnd,SampleSourceLocations,SourceBreakends,SampleInserts");

            writer.newLine();

            int lineId = 0;
            for(final LineClusterData lineData : combinedLineData)
            {
                final LineRegion primarySource = lineData.primaryRegion();
                final SvRegion combinedRegion = lineData.getCombinedPrimarySourceRegion();

                int extRegionCount = getExternalLineSampleCount(primarySource.Region);

                final String lineDefn = String.format("%d,%s,%s,%d,%d,%d,%d,%d",
                        lineId, primarySource.LineType, combinedRegion.Chromosome, combinedRegion.start(), combinedRegion.end(),
                        lineData.sampleCount(), lineData.insertRegionsCount(), extRegionCount);

                writer.write(lineDefn);
                writer.write(String.format(",%s", lineData.sampleClusterData()));
                writer.newLine();

                for(int i = 0; i < lineData.MatchedClusters.size(); ++i)
                {
                    final LineClusterData clusterData = lineData.MatchedClusters.get(i);
                    writer.write(lineDefn);
                    writer.write(String.format(",%s", clusterData.sampleClusterData()));
                    writer.newLine();
                }

                ++lineId;
            }

            writer.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write to line cluster data: {}", e.toString());
        }
    }

    private int getExternalLineSampleCount(final SvRegion region)
    {
        if(mExtLineSampleCounts.isEmpty())
            return 0;

        final Map.Entry<SvRegion,Integer> extRegion = mExtLineSampleCounts.entrySet().stream()
                .filter(x -> positionsOverlap(region.start(), region.end(),
                        x.getKey().start() - LINE_ELEMENT_PROXIMITY_DISTANCE, x.getKey().end() + LINE_ELEMENT_PROXIMITY_DISTANCE))
                .findFirst().orElse(null);

        return extRegion != null ? extRegion.getValue() : 0;
    }

    private static final String SV_DATA_FILE = "sv_data_file";
    private static final String EXT_DATA_FILE = "ext_data_file";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(SV_DATA_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(EXT_DATA_FILE, true, "External LINE data sample counts");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to write results");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        CohortLineElements cohortLineElements = new CohortLineElements(cmd);
        cohortLineElements.run();

        LOGGER.info("LINE element processing complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }




}
