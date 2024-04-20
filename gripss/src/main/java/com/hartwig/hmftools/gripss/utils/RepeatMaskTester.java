package com.hartwig.hmftools.gripss.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.gripss.GripssConfig.APP_NAME;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.gripss.rm.RepeatMaskAnnotation;
import com.hartwig.hmftools.gripss.rm.RepeatMaskAnnotator;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class RepeatMaskTester
{
    private final String mSvDataFile;
    private final RepeatMaskAnnotator mRmAnnotation;
    private final BufferedWriter mWriter;
    private final int mThreads;

    private static final String SV_DATA_FILE = "sv_data_file";
    private static final String OUTPUT_FILE = "output_file";

    public RepeatMaskTester(final ConfigBuilder configBuilder)
    {
        mRmAnnotation = new RepeatMaskAnnotator();
        mSvDataFile = configBuilder.getValue(SV_DATA_FILE);

        if(!mRmAnnotation.load(configBuilder.getValue(REPEAT_MASK_FILE), RefGenomeVersion.V37))
            System.exit(1);

        mWriter = initialiseOutput(configBuilder.getValue(OUTPUT_FILE));
        mThreads = parseThreads(configBuilder);
    }

    public void run()
    {
        List<SvData> svDataList = loadSvData();

        List<CompareTask> compareTasks = Lists.newArrayList();

        for(int i = 0; i < mThreads; ++i)
        {
            compareTasks.add(new CompareTask(i));
        }

        int taskIndex = 0;
        for(SvData svData : svDataList)
        {
            if(taskIndex >= compareTasks.size())
                taskIndex = 0;

            compareTasks.get(taskIndex).SvDataList.add(svData);

            ++taskIndex;
        }

        final List<Callable> callableList = compareTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        closeBufferedWriter(mWriter);

        GR_LOGGER.info("Gripss repeat mask comparison complete");
    }

    private List<SvData> loadSvData()
    {
        List<SvData> svDataList = Lists.newArrayList();

        try
        {
            BufferedReader fileReader = createBufferedReader(mSvDataFile);

            String line = fileReader.readLine();

            String delim = line.contains("\t") ? "\t" : ",";

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, delim);
            int sampleIndex = fieldsIndexMap.get("SampleId");
            int vcfIndex = fieldsIndexMap.get("VcfId");
            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int orientIndex = fieldsIndexMap.get("Orientation");
            int insSeqIndex = fieldsIndexMap.get("InsertSequence");
            int alignmentsIndex = fieldsIndexMap.get("InsSeqAlignments");
            int rmClassIndex = fieldsIndexMap.get("RmClass");
            int rmTypeIndex = fieldsIndexMap.get("RmType");
            int rmOrientIndex = fieldsIndexMap.get("RmOrientation");
            int rmCovIndex = fieldsIndexMap.get("RmCoverage");

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(delim, -1);

                SvData svData = new SvData(
                        values[sampleIndex], values[vcfIndex], values[chrIndex], Integer.parseInt(values[posIndex]),
                        Byte.parseByte(values[orientIndex]), values[rmClassIndex], values[rmTypeIndex], values[insSeqIndex],
                        values[alignmentsIndex].replaceAll(",", ";"),
                        Byte.parseByte(values[rmOrientIndex]), Double.parseDouble(values[rmCovIndex]));

                svDataList.add(svData);
            }

            GR_LOGGER.info("loaded {} SV data entries from file({})", svDataList.size(), mSvDataFile);
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to load SV data from file({}): {}", mSvDataFile, e.toString());
            System.exit(1);
        }

        return svDataList;
    }

    private class CompareTask implements Callable
    {
        private final int mTaskId;
        public final List<SvData> SvDataList;

        public CompareTask(int taskId)
        {
            mTaskId = taskId;
            SvDataList = Lists.newArrayList();
        }

        @Override
        public Long call()
        {
            for(int i = 0; i < SvDataList.size(); ++i)
            {
                SvData svData = SvDataList.get(i);
                processSvData(svData);

                if(i > 0 && (i % 100) == 0)
                {
                    GR_LOGGER.info("{}: processed {} SVs", mTaskId, i);
                }
            }

            GR_LOGGER.info("{}: complete", mTaskId);

            return (long)0;
        }
    }

    private void processSvData(final SvData svData)
    {
        if(svData.InsertSeqAlignments.isEmpty())
            return;

        RepeatMaskAnnotation rmAnnotation = mRmAnnotation.annotate(svData.InsertSequence, svData.InsertSeqAlignments);

        write(svData, rmAnnotation);
    }

    private BufferedWriter initialiseOutput(final String outputFilename)
    {
        try
        {
            GR_LOGGER.info("writing comparison file: {}", outputFilename);

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("SampleId,VcfId,Chromosome,Position,Orientation,Alignments,OrigRmClass,OrigRmType,OrigRmCoverage");
            writer.write(",RmId,RmRegion,RmClass,RmType,AlignmentRegion,Coverage");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void write(final SvData svData, final RepeatMaskAnnotation rmAnnotation)
    {
        try
        {
            mWriter.write(format("%s,%s,%s,%d,%d,%s,%s,%s,%.3f",
                    svData.SampleId, svData.VcfId, svData.Chromosome, svData.Position, svData.Orientation, svData.InsertSeqAlignments,
                    svData.RmClass, svData.RmType, svData.RmCoverage));

            if(rmAnnotation != null)
            {
                mWriter.write(format(",%d,%s,%s,%s,%s,%.3f",
                        rmAnnotation.RmData.Id, rmAnnotation.RmData.Region, rmAnnotation.RmData.ClassType, rmAnnotation.RmData.Repeat,
                        rmAnnotation.Alignment.toString(), rmAnnotation.Coverage));
            }
            else
            {
                mWriter.write(",-1,,,,0");
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to write annotation data: {}", e.toString());
            System.exit(1);
        }
    }

    private class SvData
    {
        // SampleId,VcfId,Chromosome,Position,Orientation,InsertSequence,RmClass,RmType,InsSeqAlignments,RmOrientation,RmCoverage
        public final String SampleId;
        public final String VcfId;
        public final String Chromosome;
        public final int Position;
        public final byte Orientation;
        public final String RmClass;
        public final String RmType;
        public final String InsertSequence;
        public final String InsertSeqAlignments;
        public final byte RmOrientation;
        public final double RmCoverage;

        public SvData(
                final String sampleId, final String vcfId, final String chromosome, final int position, final byte orientation,
                final String rmClass, final String rmType, final String insSeq, final String alignments, final byte rmOrientation,
                final double rmCoverage)
        {
            SampleId = sampleId;
            VcfId = vcfId;
            Chromosome = chromosome;
            Position = position;
            Orientation = orientation;
            RmClass = rmClass;
            RmType = rmType;
            InsertSequence = insSeq;
            InsertSeqAlignments = alignments;
            RmOrientation = rmOrientation;
            RmCoverage = rmCoverage;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addPath(SV_DATA_FILE, true, "File with SV data");
        configBuilder.addPath(OUTPUT_FILE, true, "Output matching file");
        addThreadOptions(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GR_LOGGER.info("running Gripss repeat-mask annotation tester");

        RepeatMaskTester repeatMaskTester = new RepeatMaskTester(configBuilder);
        repeatMaskTester.run();
    }
}
