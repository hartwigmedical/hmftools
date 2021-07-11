package com.hartwig.hmftools.svtools.cohort;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SineBreakendFinder
{
    private static final Logger LOGGER = LogManager.getLogger(SineBreakendFinder.class);

    private final List<String> mSampleIds;
    private final Map<String,List<MatchedSineElement>> mMatchedSineElements;

    private final DatabaseAccess mDbAccess;

    private final List<String> mRmTypes;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String REPEAT_MASKER_FILE = "repeat_masker_file";
    private static final String RM_TYPES = "rm_types";

    private static final int SINE_ELEMENT_DISTANCE = 1000;
    private static final int REPEAT_SUB_LENGTH = 4;

    private final String mOutputDir;

    public SineBreakendFinder(final CommandLine cmd)
    {
        String sampleIdsFile = cmd.getOptionValue(SAMPLE_ID_FILE);
        mSampleIds = ConfigUtils.loadSampleIdFile(sampleIdsFile);
        LOGGER.info("loaded {} samples from file()", mSampleIds.size(), cmd.getOptionValue(SAMPLE_ID_FILE));

        mRmTypes = Lists.newArrayList();

        if(cmd.hasOption(RM_TYPES))
        {
            mRmTypes.addAll(Arrays.stream(cmd.getOptionValue(RM_TYPES).split(";", -1)).map(x -> x).collect(Collectors.toList()));
        }

        LOGGER.info("loaded {} repeat masker types", mRmTypes.size());

        mMatchedSineElements = Maps.newHashMap();
        loadRepeatMaskerDataFile(cmd.getOptionValue(REPEAT_MASKER_FILE));

        mOutputDir = parseOutputDir(cmd);
        mDbAccess = createDatabaseAccess(cmd);
    }

    public void run()
    {
        if(mMatchedSineElements.isEmpty())
        {
            LOGGER.info("no matched SINE elements found");
            return;
        }

        if(mDbAccess == null)
        {
            LOGGER.error("no database connection");
            return;
        }

        if(mSampleIds.isEmpty())
        {
            mSampleIds.addAll(mDbAccess.readPurpleSampleListPassingQC(MIN_SAMPLE_PURITY));
            LOGGER.info("loaded {} samples from database", mSampleIds.size());
        }

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            final List<StructuralVariantData> svRecords = mDbAccess.readStructuralVariantData(sampleId);

            if(svRecords.isEmpty())
                continue;

            svRecords.forEach(x -> checkVariant(sampleId, x));

            if(i > 0 && (i % 10) == 0)
            {
                int matchedBreakends = mMatchedSineElements.values().stream().mapToInt(x -> x.stream().mapToInt(y -> y.Breakends.size()).sum()).sum();
                LOGGER.info("processed {} samples, matched breakends({})", i, matchedBreakends);
            }
        }

        writeResults();
    }

    private void checkVariant(final String sampleId, final StructuralVariantData svData)
    {
        if(!(svData.filter().isEmpty() || svData.filter().equals(PASS) || svData.filter().equals(INFERRED)))
            return;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && (svData.type() == SGL || svData.type() == INF))
                break;

            final String chromosome = se == SE_START ? svData.startChromosome() : svData.endChromosome();
            final byte orient = se == SE_START ? svData.startOrientation() : svData.endOrientation();
            final int position = se == SE_START ? svData.startPosition() : svData.endPosition();

            final List<MatchedSineElement> rmElements = mMatchedSineElements.get(chromosome);

            if(rmElements == null)
                continue;

            for(MatchedSineElement rmElement : rmElements)
            {
                if(rmElement.isPositionWithin(position))
                {
                    rmElement.Breakends.add(new SvBreakendData(sampleId, svData.id(), position, se == SE_START, orient, svData.type()));
                }
            }
        }
    }

    private void loadRepeatMaskerDataFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            int rmIdIndex = fieldsIndexMap.get("RmId");
            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("PosStart");
            int posEndIndex = fieldsIndexMap.get("PosEnd");
            int strandIndex = fieldsIndexMap.get("Strand");
            int repeatIndex = fieldsIndexMap.get("MatchingRepeat");
            int classIndex = fieldsIndexMap.get("ClassFamily");

            Map<String,List<RepeatMaskerData>> rmDataSubTypeMap = Maps.newHashMap();

            for(final String line : fileContents)
            {
                final String[] items = line.split(",");

                final String matchingRepeat = items[repeatIndex];

                if(!mRmTypes.isEmpty() && !mRmTypes.contains(matchingRepeat))
                    continue;

                int rmId = Integer.parseInt(items[rmIdIndex]);

                final String chromosome = items[chrIndex].replace("chr", "");

                if(!HumanChromosome.contains(chromosome))
                    continue;

                final int[] positions =
                        new int[] { Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]) };

                byte strand = items[strandIndex].equals("+") ? POS_ORIENT : NEG_ORIENT;

                RepeatMaskerData rmData = new RepeatMaskerData(
                        rmId, new BaseRegion(chromosome, positions), strand, items[classIndex], matchingRepeat);

                final String subRepeat = rmData.subRepeat(REPEAT_SUB_LENGTH);

                List<RepeatMaskerData> rmDataList = rmDataSubTypeMap.get(subRepeat);

                if(rmDataList == null)
                {
                    rmDataSubTypeMap.put(subRepeat, Lists.newArrayList(rmData));
                }
                else
                {
                    rmDataList.add(rmData);
                }
            }

            for(List<RepeatMaskerData> rmDataList : rmDataSubTypeMap.values())
            {
                RepeatMaskerData lastRmData = null;

                for(RepeatMaskerData rmData : rmDataList)
                {
                    if(lastRmData != null && areMatchingPair(lastRmData, rmData))
                    {
                        MatchedSineElement matchedRepeat = new MatchedSineElement(lastRmData, rmData);

                        LOGGER.debug("matching RM elements: first({}) second({})", lastRmData, rmData);

                        final String chromosome = rmData.Region.Chromosome;

                        List<MatchedSineElement> matchedElements = mMatchedSineElements.get(chromosome);

                        if(matchedElements == null)
                        {
                            matchedElements = Lists.newArrayList(matchedRepeat);
                            mMatchedSineElements.put(chromosome, matchedElements);
                        }
                        else
                        {
                            matchedElements.add(matchedRepeat);
                        }
                    }

                    lastRmData = rmData;
                }
            }

            LNX_LOGGER.info("loaded {} repeat-masker items, pairs({}) items from file: {}",
                    fileContents.size(), mMatchedSineElements.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read repeat-masker CSV file({})", filename);
        }
    }

    private boolean areMatchingPair(final RepeatMaskerData first, final RepeatMaskerData second)
    {
        if(!first.Region.Chromosome.equals(second.Region.Chromosome))
            return false;

        if(first.Strand == second.Strand)
            return false;

        if(first.Region.end() >= second.Region.start())
            return false;

        if(!first.isProximateElement(second, SINE_ELEMENT_DISTANCE))
            return false;

        return true;
    }

    private void writeResults()
    {
        try
        {
            String outputFileName = mOutputDir + "SINE_ELEMENT_BREAKENDS.csv";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId,Chromosome,RmIds");
            writer.write(",FirstRepeat,FirstOrient,FirstPosStart,FirstPosEnd,SecondRepeat,SecondOrient,SecondPosStart,SecondPosEnd");
            writer.write(",SvId,Type,IsStart,Position,Orientation");

            writer.newLine();

            for(Map.Entry<String,List<MatchedSineElement>> entry : mMatchedSineElements.entrySet())
            {
                final String chromosome = entry.getKey();

                for(MatchedSineElement rmElement : entry.getValue())
                {
                    for(SvBreakendData breakend : rmElement.Breakends)
                    {
                        writer.write(String.format("%s,%s,%s", breakend.SampleId, chromosome, rmElement.combinedRmId()));

                        writer.write(String.format(",%s,%d,%d,%d,%s,%d,%d,%d",
                                rmElement.RmStart.Repeat, rmElement.RmStart.Strand, rmElement.RmStart.Region.start(), rmElement.RmStart.Region.end(),
                                rmElement.RmEnd.Repeat, rmElement.RmEnd.Strand, rmElement.RmEnd.Region.start(), rmElement.RmEnd.Region.end()));

                        writer.write(String.format(",%d,%s,%s,%d,%d",
                                breakend.SvId, breakend.Type, breakend.IsStart, breakend.Position, breakend.Orientation));

                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write to line cluster data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(SAMPLE_ID_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(RM_TYPES, true, "External LINE data sample counts");
        options.addOption(REPEAT_MASKER_FILE, true, "Polymorphic LINE data file");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        DatabaseAccess.addDatabaseCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SineBreakendFinder cohortLineElements = new SineBreakendFinder(cmd);
        cohortLineElements.run();

        LOGGER.info("SINE breakend matching complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
