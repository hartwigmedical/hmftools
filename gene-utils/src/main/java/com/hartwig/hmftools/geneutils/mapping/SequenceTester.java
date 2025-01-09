package com.hartwig.hmftools.geneutils.mapping;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class SequenceTester
{
    private final SeqTestConfig mConfig;
    private final BlastNSeqTester mBlastnTester;

    private final BwaSeqTester mBwaSeqTester;

    protected static final String FLD_SEQ_ID = "SeqId";
    protected static final String FLD_SEQUENCE = "Sequence";

    public SequenceTester(final ConfigBuilder configBuilder)
    {
        mConfig = new SeqTestConfig(configBuilder);

        mBlastnTester = new BlastNSeqTester(mConfig, configBuilder);
        mBwaSeqTester = new BwaSeqTester(mConfig, configBuilder);
    }

    public void run()
    {
        if(mConfig.InputFile == null || !Files.exists(Paths.get(mConfig.InputFile)))
        {
            GU_LOGGER.error("invalid input file");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        List<SequenceInfo> sequenceInfos = loadSequenceInfo();

        mBwaSeqTester.run(sequenceInfos);

        mBlastnTester.run(sequenceInfos);

        GU_LOGGER.info("sequence tester complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<SequenceInfo> loadSequenceInfo()
    {
        List<SequenceInfo> sequenceInfos = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(mConfig.InputFile));

            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            Integer chrIndex = fieldsIndexMap.getOrDefault(FLD_CHROMOSOME, null);
            Integer posStartIndex = fieldsIndexMap.getOrDefault(FLD_POSITION_START, null);
            Integer posEndIndex = fieldsIndexMap.getOrDefault(FLD_POSITION_END, null);
            Integer seqIndex = fieldsIndexMap.getOrDefault(FLD_SEQUENCE, null);

            int seqId = 0;
            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM);

                String chromosome = chrIndex != null ? values[chrIndex] : "";
                String startString = posStartIndex != null ? values[posStartIndex] : "";
                String endString = posEndIndex != null ? values[posEndIndex] : "";
                String sequence = seqIndex != null ? values[seqIndex] : "";

                ChrBaseRegion region = null;
                if(chromosome.length() > 0 && startString.length() > 0 && endString.length() > 0)
                {
                    region = new ChrBaseRegion(chromosome, Integer.parseInt(startString), Integer.parseInt(endString));
                }

                if(sequence == null || sequence.length() == 0)
                {
                    sequence = mConfig.RefGenome.getBaseString(chromosome, region.start(), region.end());
                }

                sequenceInfos.add(new SequenceInfo(seqId, sequence, region));
                seqId++;
            }

            GU_LOGGER.info("loaded {} sequences from file({})", sequenceInfos.size(), mConfig.InputFile);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error reading input file({}): {}", mConfig.InputFile, e.toString());
            System.exit(1);
        }

        return sequenceInfos;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        SeqTestConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SequenceTester sequenceTester = new SequenceTester(configBuilder);

        sequenceTester.run();
    }
}
