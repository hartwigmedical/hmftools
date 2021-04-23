package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.REF_GENOME_FILE;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.LINE_ELEMENT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.fromString;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class CohortLineElements
{
    private static final Logger LOGGER = LogManager.getLogger(CohortLineElements.class);

    private final Map<String,Map<Integer,LineClusterData>> mSampleClusterLineData;
    private final Map<BaseRegion,Integer> mExtLineSampleCounts;
    private final Map<String,List<RepeatMaskerData>> mChrRepeatMaskerData;
    private final List<BaseRegion> mKnownLineElements;
    private final List<BaseRegion> mPolymorphicLineElements;
    private IndexedFastaSequenceFile mRefGenomeFile;

    private final String mOutputDir;
    private final String mSvDataFile;
    private final String mExtDataFile;
    private final String mRepeatMaskerDataFile;
    private final String mPolymorphicDataFile;
    private final String mKnownLineElementsFile;
    public final boolean mWriteRefGenomeLineBases;

    public CohortLineElements(final CommandLine cmd)
    {
        mSampleClusterLineData = Maps.newHashMap();
        mExtLineSampleCounts = Maps.newHashMap();
        mChrRepeatMaskerData = Maps.newHashMap();
        mKnownLineElements = Lists.newArrayList();
        mPolymorphicLineElements = Lists.newArrayList();

        mOutputDir = parseOutputDir(cmd);
        mSvDataFile = cmd.getOptionValue(SV_DATA_FILE);
        mExtDataFile = cmd.getOptionValue(EXT_DATA_FILE);
        mPolymorphicDataFile = cmd.getOptionValue(POLYMORPHIC_DATA_FILE);
        mRepeatMaskerDataFile = cmd.getOptionValue(REPEAT_MASKER_DATA_FILE);
        mKnownLineElementsFile = cmd.getOptionValue(KNOWN_DATA_FILE);
        mWriteRefGenomeLineBases = cmd.hasOption(WRITE_LINE_SEQUENCES);
        mRefGenomeFile = null;

        try
        {
            if(cmd.hasOption(REF_GENOME_FILE))
            {
                mRefGenomeFile = new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME_FILE)));
            }
        }
        catch (Exception e)
        {
            LNX_LOGGER.error("failed to load ref genome: {}", e.toString());
        }
    }

    public void run()
    {
        if(mSvDataFile == null)
            return;

        loadLineElementsFile(mSvDataFile);
        loadExternalLineDataFile(mExtDataFile);
        loadRepeatMaskerLineDataFile(mRepeatMaskerDataFile);
        loadKnownLineElements(mKnownLineElementsFile);
        loadPolymorphicLineElements(mPolymorphicDataFile);

        produceResults();

        if(mWriteRefGenomeLineBases && !mChrRepeatMaskerData.isEmpty())
        {
            writeRefGeneLineBases();
        }
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

    public void loadKnownLineElements(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileContents)
            {
                // parse CSV data
                String[] items = line.split(",");

                final BaseRegion lineRegion = new BaseRegion(
                        RG_VERSION.versionedChromosome(items[fieldsIndexMap.get("Chromosome")]),
                        Integer.parseInt(items[fieldsIndexMap.get("PosStart")]),
                        Integer.parseInt(items[fieldsIndexMap.get("PosEnd")]));

                mKnownLineElements.add(lineRegion);
            }

            LNX_LOGGER.info("loaded {} known line elements from file: {}", mKnownLineElements.size(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    public void loadPolymorphicLineElements(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            final String header = fileContents.get(0);
            fileContents.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileContents)
            {
                // parse CSV data
                String[] items = line.split(",");

                int pos1 = Integer.parseInt(items[fieldsIndexMap.get("PositivePosition")]);
                int pos2 = Integer.parseInt(items[fieldsIndexMap.get("NegativePosition")]);

                final BaseRegion lineRegion = new BaseRegion(
                        RG_VERSION.versionedChromosome(items[fieldsIndexMap.get("Chromosome")]),
                        Math.min(pos1, pos2), Math.max(pos1, pos2));

                mPolymorphicLineElements.add(lineRegion);
            }

            LNX_LOGGER.info("loaded {} polymorphic line elements from file: {}", mPolymorphicLineElements.size(), filename);
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

                mExtLineSampleCounts.put(new BaseRegion(chromosome, positions), sampleCount);
            }

            LNX_LOGGER.info("loaded {} external line data items from file: {}", mExtLineSampleCounts.size(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    private void loadRepeatMaskerLineDataFile(final String filename)
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

            String currentChr = "";
            List<RepeatMaskerData> rmDataList = null;
            int itemCount = 0;

            for(final String line : fileContents)
            {
                final String[] items = line.split(",");

                double rmIdDbl = Double.parseDouble(items[rmIdIndex]);
                final int rmId = (int)rmIdDbl;

                final String chromosome = items[chrIndex];

                final int[] positions =
                        new int[] { Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]) };

                final byte strand = items[strandIndex].equals("+") ? POS_ORIENT : NEG_ORIENT;

                RepeatMaskerData rmData = new RepeatMaskerData(rmId, new BaseRegion(chromosome, positions), strand);
                ++itemCount;

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;

                    rmDataList = mChrRepeatMaskerData.get(chromosome);

                    if(rmDataList == null)
                    {
                        rmDataList = Lists.newArrayList();
                        mChrRepeatMaskerData.put(chromosome, rmDataList);
                    }
                }

                rmDataList.add(rmData);
            }

            LNX_LOGGER.info("loaded {} repeat-masker line data items from file: {}", itemCount, filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read line element CSV file({})", filename);
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
                            sampleId, clusterId, new BaseRegion(chromosomes[se], positions[se], positions[se]), lineTypes[se]);
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

            writer.write("LineId,Type,Chromosome,PosStart,PosEnd,SampleCount,TotalInserts,PcawgSampleCount");
            writer.write(",KnownPosStart,KnownPosEnd,RmId,RmPosStart,RmPosEnd,RmStrand,PolymorphPosStart,PolymorphPosEnd");
            writer.write(",SampleId,ClusterId,SamplePosStart,SamplePosEnd,SampleSourceLocations,SourceBreakends,SampleInserts");

            writer.newLine();

            int lineId = 0;
            for(final LineClusterData lineData : combinedLineData)
            {
                final LineRegion primarySource = lineData.primaryRegion();
                final BaseRegion combinedRegion = lineData.getCombinedPrimarySourceRegion();

                final BaseRegion knownLineRegion = findKnownLineRegion(primarySource.Region);

                final BaseRegion pmLineRegion = knownLineRegion == null ? findPolymorphicLineRegion(primarySource.Region) : null;

                RepeatMaskerData rmData = findRepeatMaskerMatch(primarySource.Region, knownLineRegion);

                final String rmDataStr = rmData != null ?
                        String.format("%d,%d,%d,%d", rmData.RmId, rmData.Region.start(), rmData.Region.end(), rmData.Strand) : "-1,-1,-1,0";

                int extRegionCount = getExternalLineSampleCount(primarySource.Region);

                final String lineDefn = String.format("%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%s,%d,%d",
                        lineId, primarySource.LineType, combinedRegion.Chromosome, combinedRegion.start(), combinedRegion.end(),
                        lineData.sampleCount(), lineData.insertRegionsCount(), extRegionCount,
                        knownLineRegion != null ? knownLineRegion.start() : -1, knownLineRegion != null ? knownLineRegion.end() : -1,
                        rmDataStr, pmLineRegion != null ? pmLineRegion.start() : -1, pmLineRegion != null ? pmLineRegion.end() : -1);

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

    private BaseRegion findKnownLineRegion(final BaseRegion lineRegion)
    {
        for(BaseRegion knownRegion : mKnownLineElements)
        {
            if(!lineRegion.Chromosome.equals(knownRegion.Chromosome))
                continue;

            int[] knownRegionRange = { knownRegion.start() - LINE_ELEMENT_PROXIMITY_DISTANCE, knownRegion.end() + LINE_ELEMENT_PROXIMITY_DISTANCE };

            if(positionsOverlap(lineRegion.start(), lineRegion.end(), knownRegionRange[SE_START], knownRegionRange[SE_END]))
            {
                return knownRegion;
            }
        }

        return null;
    }

    private static int regionMidpoint(final BaseRegion region) { return (region.end() + region.end())  / 2; }

    private BaseRegion findPolymorphicLineRegion(final BaseRegion lineRegion)
    {
        BaseRegion closestPmRegion = null;
        int closestDistance = -1;

        for(BaseRegion pmRegion : mPolymorphicLineElements)
        {
            if(!lineRegion.Chromosome.equals(pmRegion.Chromosome))
                continue;

            int[] pmRegionRange = { pmRegion.start() - LINE_ELEMENT_PROXIMITY_DISTANCE, pmRegion.end() + LINE_ELEMENT_PROXIMITY_DISTANCE };

            if(positionsOverlap(lineRegion.start(), lineRegion.end(), pmRegionRange[SE_START], pmRegionRange[SE_END]))
            {
                int distance = abs(regionMidpoint(lineRegion) - regionMidpoint(pmRegion));

                if(closestPmRegion == null || distance < closestDistance)
                {
                    closestPmRegion = pmRegion;
                    closestDistance = distance;
                }
            }
        }

        return closestPmRegion;
    }

    private static final int INTACT_LINE_ELEMENT_LENGTH = 5000;

    private RepeatMaskerData findRepeatMaskerMatch(final BaseRegion lineRegion, final BaseRegion knownLineRegion)
    {
        final List<RepeatMaskerData> rmDataList = mChrRepeatMaskerData.get(lineRegion.Chromosome);

        if(rmDataList == null)
            return null;

        RepeatMaskerData closestRmData = null;
        int closestDistance = -1;

        for(RepeatMaskerData rmData : rmDataList)
        {
            if(knownLineRegion != null && knownLineRegion.overlaps(rmData.Region))
            {
                return rmData;
            }

            if(rmData.Region.length() < INTACT_LINE_ELEMENT_LENGTH)
                continue;

            // teh line element is expected to be up-stream of the SV activity
            final int[] rmMatchLimits = {0, 0};

            if(rmData.Strand == POS_ORIENT)
            {
                rmMatchLimits[SE_START] = rmData.Region.start() - LINE_ELEMENT_PROXIMITY_DISTANCE;
                rmMatchLimits[SE_END] = rmData.Region.end() + LINE_ELEMENT_PROXIMITY_DISTANCE * 4;
            }
            else
            {
                rmMatchLimits[SE_START] = rmData.Region.start() - LINE_ELEMENT_PROXIMITY_DISTANCE * 4;
                rmMatchLimits[SE_END] = rmData.Region.end() + LINE_ELEMENT_PROXIMITY_DISTANCE;
            }

            if(positionsOverlap(lineRegion.start(), lineRegion.end(), rmMatchLimits[SE_START], rmMatchLimits[SE_END]))
            {
                int distance = abs(regionMidpoint(lineRegion) - regionMidpoint(rmData.Region));
                if(closestRmData == null || distance < closestDistance)
                {
                    closestRmData = rmData;
                    closestDistance = distance;
                }
            }
        }

        return closestRmData;
    }

    private int getExternalLineSampleCount(final BaseRegion region)
    {
        if(mExtLineSampleCounts.isEmpty())
            return 0;

        final Map.Entry<BaseRegion,Integer> extRegion = mExtLineSampleCounts.entrySet().stream()
                .filter(x -> positionsOverlap(region.start(), region.end(),
                        x.getKey().start() - LINE_ELEMENT_PROXIMITY_DISTANCE, x.getKey().end() + LINE_ELEMENT_PROXIMITY_DISTANCE))
                .findFirst().orElse(null);

        return extRegion != null ? extRegion.getValue() : 0;
    }

    private void writeRefGeneLineBases()
    {
        if(mChrRepeatMaskerData.isEmpty() || mRefGenomeFile == null)
            return;

        try
        {
            String outputFileName = mOutputDir + "line_ref_bases.fasta";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            for(Map.Entry<String,List<RepeatMaskerData>> entry : mChrRepeatMaskerData.entrySet())
            {
                final String chromosome = entry.getKey();

                for(final RepeatMaskerData rmData : entry.getValue())
                {
                    if(rmData.Region.length() < 5000)
                        continue;

                    int[] coords = new int[] { rmData.Region.start(), rmData.Region.end() };

                    if(rmData.Strand == POS_ORIENT)
                        coords[SE_END] += 5000;
                    else
                        coords[SE_START] -= 5000;

                    final String refBases = mRefGenomeFile.getSubsequenceAt(chromosome, rmData.Region.start(), rmData.Region.end()).getBaseString();

                    if(refBases == null || refBases.isEmpty())
                        continue;

                    writer.write(String.format(">%d: %s:%d-%d", rmData.RmId, chromosome, coords[SE_START], coords[SE_END]));
                    writer.newLine();

                    int index = 0;
                    while(index < refBases.length())
                    {
                        writer.write(refBases.substring(index, min(index + 80, refBases.length())));
                        writer.newLine();
                        index += 80;
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write to line ref-genome bases: {}", e.toString());
        }
    }

    private static final String SV_DATA_FILE = "sv_data_file";
    private static final String EXT_DATA_FILE = "ext_data_file";
    private static final String POLYMORPHIC_DATA_FILE = "polymorphic_data_file";
    private static final String KNOWN_DATA_FILE = "known_line_elements_file";
    private static final String REPEAT_MASKER_DATA_FILE = "repeat_masker_data_file";
    private static final String WRITE_LINE_SEQUENCES = "write_line_seq";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(SV_DATA_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(EXT_DATA_FILE, true, "External LINE data sample counts");
        options.addOption(POLYMORPHIC_DATA_FILE, true, "Polymorphic LINE data file");
        options.addOption(KNOWN_DATA_FILE, true, "Known LINE elements file");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(REPEAT_MASKER_DATA_FILE, true, "Path to repeat masker data for LINE elements");
        options.addOption(WRITE_LINE_SEQUENCES, false, "Write ref genome LINE element sequences");
        options.addOption(REF_GENOME_FILE, true, "Path to the indexed ref genome fasta file");
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
