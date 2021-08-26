package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.log;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BinderConfig;
import com.hartwig.hmftools.neo.bind.BlosumMapping;
import com.hartwig.hmftools.neo.bind.RandomPeptideConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

// routine for finding the most similar peptide from a collection using Blosum amino-acid correlations
public class PeptideSimilarities
{
    private final String mPeptidesFile;

    private final BlosumMapping mBlosumMapping;
    private BufferedWriter mWriter;

    private static final String PEPTIDES_FILE = "peptides_file";

    public PeptideSimilarities(final CommandLine cmd)
    {
        mPeptidesFile = cmd.getOptionValue(PEPTIDES_FILE);
        mBlosumMapping = new BlosumMapping();

        String outputDir = parseOutputDir(cmd);
        String outputFilename = BinderConfig.formFilename("peptide_similarities", outputDir, cmd.getOptionValue(OUTPUT_ID));
        mWriter = initialiseWriter(outputFilename);
    }

    public void run()
    {
        NE_LOGGER.info("loading MCF random peptide predictions from file({})", mPeptidesFile);

        final Map<String,Map<Integer,List<BindData>>> allelePeptideData = Maps.newHashMap();

        if(!loadBindData(mPeptidesFile, true, Lists.newArrayList(), Lists.newArrayList(), allelePeptideData))
            return;

        int peptideCount = allelePeptideData.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
        NE_LOGGER.info("finding top similarities for {} peptides", peptideCount);

        Set<String> processedPeptides = Sets.newHashSet();

        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : allelePeptideData.entrySet())
        {
            final Map<Integer, List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            NE_LOGGER.info("processing allele({}) with {} peptides",
                    alleleEntry.getKey(), pepLenBindDataMap.values().stream().mapToInt(x -> x.size()).sum());

            for(Map.Entry<Integer, List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<BindData> bindDataList = pepLenEntry.getValue();

                for(BindData bindData : bindDataList)
                {
                    if(processedPeptides.contains(bindData.Peptide))
                        continue;

                    processedPeptides.add(bindData.Peptide);

                    findPeptideSimilarity(bindData, allelePeptideData);
                }
            }
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("similarity search complete");
    }

    public static BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Peptide,TopSimilarity,TopPeptide");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
            return null;
        }
    }

    private void findPeptideSimilarity(final BindData bindData, final Map<String,Map<Integer,List<BindData>>> allelePeptideData)
    {
        String topPeptide = "";
        double topSimiliarity = 0;

        double selfSimilarity = mBlosumMapping.calcSequenceBlosumScore(bindData.Peptide);

        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : allelePeptideData.entrySet())
        {
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            List<BindData> otherBindDataList = pepLenBindDataMap.get(bindData.peptideLength());

            if(otherBindDataList == null)
                continue;

            for(BindData otherBindData : otherBindDataList)
            {
                if(otherBindData.Peptide.equals(bindData.Peptide))
                    continue;

                double similarity = mBlosumMapping.calcSequenceBlosumScore(otherBindData.Peptide, bindData.Peptide);

                if(similarity > topSimiliarity)
                {
                    topPeptide = otherBindData.Peptide;
                    topSimiliarity = similarity;
                }
            }
        }

        try
        {
            mWriter.write(String.format("%s,%s,%6.3e,%s",
                    bindData.Allele, bindData.Peptide, log(topSimiliarity/selfSimilarity), topPeptide));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write similarity data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(PEPTIDES_FILE, true, "MCF predictions file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PeptideSimilarities peptideSimilarities = new PeptideSimilarities(cmd);
        peptideSimilarities.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
