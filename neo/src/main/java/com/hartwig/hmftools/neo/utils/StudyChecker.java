package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.stats.FisherExactTest;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

// looks for outlier studies in allele-peptide binding training data
public class StudyChecker
{
    private final Map<String,List<BindStudyData>> mAlleleDataMap;
    private final Set<Integer> mPositions;
    private final Set<String> mStudies;
    private final String mOutputDir;
    private final BufferedWriter mWriter;

    private final FisherExactTest mFisherEt;

    private static final String STUDY_DATA_FILE = "study_data_file";

    public StudyChecker(final CommandLine cmd)
    {
        String dataFile = cmd.getOptionValue(STUDY_DATA_FILE);

        mOutputDir = parseOutputDir(cmd);

        mWriter = initProbWriter();

        mAlleleDataMap = Maps.newHashMap();
        mPositions = Sets.newHashSet();
        mStudies = Sets.newHashSet();

        loadTrainingData(dataFile);

        mFisherEt = new FisherExactTest();
        mFisherEt.initialise(1000000);
    }

    public void run()
    {
        int maxPosition = mPositions.stream().mapToInt(x -> x.intValue()).max().orElse(0);

        for(Map.Entry<String,List<BindStudyData>> entry : mAlleleDataMap.entrySet())
        {
            String allele = entry.getKey();

            NE_LOGGER.info("evaluating allele({})", allele);

            for(int pos = 0; pos <= maxPosition; ++pos)
            {
                evaluate(allele, entry.getValue(), pos);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private BufferedWriter initProbWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputDir + "study_check_probs.csv", false);

            writer.write("Allele,Position,AminoAcid,Study,Expected,Both,Prob");
            writer.write(",TotalBinds,WithStudy,WithAminoAcid,WithAminoAcidNoStudy,WithStudyNoAcid,Neither");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write study prob data: {}", e.toString());
        }

        return null;
    }

    private static final double LOW_PROB = 1e-6;

    private void evaluate(final String allele, final List<BindStudyData> studyDataList, int position)
    {
        try
        {
            for(Character aminoAcid : AMINO_ACIDS)
            {
                for(String study : mStudies)
                {
                    int totalBinds = 0;
                    int withStudy = 0;
                    int withAminoAcid = 0;
                    int withAminoAcidNoStudy = 0;
                    int withStudyNoAcid = 0;
                    int both = 0;
                    int neither = 0;

                    for(BindStudyData data : studyDataList)
                    {
                        if(data.Position != position)
                            continue;

                        totalBinds += data.BindCount;

                        if(data.AminoAcid == aminoAcid && data.Study.equals(study))
                        {
                            withAminoAcid += data.BindCount;
                            withStudy += data.BindCount;
                            both += data.BindCount;
                        }
                        else if(data.AminoAcid == aminoAcid && !data.Study.equals(study))
                        {
                            withAminoAcid += data.BindCount;
                            withAminoAcidNoStudy += data.BindCount;
                        }
                        else if(data.AminoAcid != aminoAcid && data.Study.equals(study))
                        {
                            withStudyNoAcid += data.BindCount;
                            withStudy += data.BindCount;
                        }
                        else
                        {
                            neither += data.BindCount;
                        }
                    }

                    double expected = (withAminoAcid * withStudy) / (double) totalBinds;

                    double prob = mFisherEt.calc(both, withAminoAcidNoStudy, withStudyNoAcid, neither, expected);

                    if(prob < LOW_PROB)
                    {
                        NE_LOGGER.debug(String.format("allele(%s) position(%d) aminoAcids(%s) study(%s) expected(%.3f) prob(%4.3e)",
                                allele, position, aminoAcid, study, expected, prob));

                        mWriter.write(String.format("%s,%s,%s,%s,%.2f,%d,%4.3e,%d,%d,%d,%d,%d,%d",
                                allele, position, aminoAcid, study, expected, both, prob, totalBinds,
                                withStudy, withAminoAcid, withAminoAcidNoStudy, withStudyNoAcid, neither));
                        mWriter.newLine();
                    }
                }
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write study prob data: {}", e.toString());
        }
    }

    private void loadTrainingData(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int studyIndex = fieldsIndexMap.get("Study");
            int posIndex = fieldsIndexMap.get("Position");
            int aaIndex = fieldsIndexMap.get("AminoAcid");
            int bcIndex = fieldsIndexMap.get("BindCount");

            String currentAllele = "";
            List<BindStudyData> currentBindList = null;

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = items[alleleIndex];

                if(!allele.equals(currentAllele))
                {
                    currentAllele = allele;
                    currentBindList = Lists.newArrayList();
                    mAlleleDataMap.put(allele, currentBindList);
                }

                BindStudyData data = new BindStudyData(
                        allele, items[studyIndex], items[aaIndex].charAt(0), Integer.parseInt(items[posIndex]), Integer.parseInt(items[bcIndex]));

                mPositions.add(data.Position);
                mStudies.add(data.Study);

                currentBindList.add(data);
            }

            NE_LOGGER.info("loaded {} alleles, {} positions, {} studies with {} study data items from file({})",
                    mAlleleDataMap.size(), mPositions.size(), mStudies.size(),
                    mAlleleDataMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read study data file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(STUDY_DATA_FILE, true, "Study binding data");
        addOutputDir(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        StudyChecker studyChecker = new StudyChecker(cmd);
        studyChecker.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private class BindStudyData
    {
        public final String Allele;
        public final String Study;
        public final Character AminoAcid;
        public final int Position;
        public final int BindCount;

        public BindStudyData(final String allele, final String study, final char aminoAcid, final int position, final int bindCount)
        {
            Allele = allele;
            Study = study;
            AminoAcid = aminoAcid;
            Position = position;
            BindCount = bindCount;
        }
    }
}
