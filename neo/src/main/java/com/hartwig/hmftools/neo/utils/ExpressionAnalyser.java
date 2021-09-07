package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_LIKE_RANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.TranscriptExpression;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class ExpressionAnalyser
{
    private final TranscriptExpression mTranscriptExpression;
    private final Map<String,PeptideSearchData> mPeptideSearchDataMap;
    private final List<PeptideData> mProteomePeptides;
    private final List<PeptideData> mBinderPeptides;

    private final Set<String> mValidationAlleles;

    private final BufferedWriter mWriter;

    private static final String PROTEOME_PEPTIDES_FILE = "proteome_binders_file";
    private static final String VALIDATION_PEPTIDES_FILE = "validation_file";
    private static final String PEPTIDE_SEARCH_FILE = "peptide_search_file";

    public ExpressionAnalyser(final CommandLine cmd)
    {
        mTranscriptExpression = new TranscriptExpression(cmd.getOptionValue(IMMUNE_EXPRESSION_FILE));

        mBinderPeptides = loadPeptideData(cmd.getOptionValue(VALIDATION_PEPTIDES_FILE), "Validation");
        mValidationAlleles = mBinderPeptides.stream().map(x -> x.Allele).collect(Collectors.toSet());

        mProteomePeptides = loadPeptideData(cmd.getOptionValue(PROTEOME_PEPTIDES_FILE), "Proteome");

        mPeptideSearchDataMap = loadPeptideSearchData(cmd.getOptionValue(PEPTIDE_SEARCH_FILE));

        String outputDir = parseOutputDir(cmd);
        String outputFile = outputDir + "peptide_expression.csv";
        mWriter = initialiseWriter(outputFile);
    }

    public void run()
    {
        if(!mTranscriptExpression.hasValidData() || mPeptideSearchDataMap.isEmpty() || mProteomePeptides.isEmpty() || mBinderPeptides.isEmpty())
        {
            NE_LOGGER.error("failed to initialise");
            System.exit(1);
        }

        NE_LOGGER.info("finding expression data for {} validation peptides and {} proteome peptides",
                mBinderPeptides.size(), mProteomePeptides.size());

        for(PeptideData peptideData : mBinderPeptides)
        {
            // find transcripts based on search results
            PeptideSearchData searchData = mPeptideSearchDataMap.get(peptideData.Peptide);
            if(searchData == null)
            {
                NE_LOGGER.warn("peptide({}) missing search data");
                continue;
            }

            peptideData.Transcripts.addAll(searchData.Transcripts);

            for(String transName : peptideData.Transcripts)
            {
                peptideData.addTpm(mTranscriptExpression.getExpression(transName));
            }

            writePeptideData(peptideData);
        }

        for(PeptideData peptideData : mProteomePeptides)
        {
            if(!mValidationAlleles.contains(peptideData.Allele))
                continue;

            for(String transName : peptideData.Transcripts)
            {
                peptideData.addTpm(mTranscriptExpression.getExpression(transName));
            }

            writePeptideData(peptideData);
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("peptide expression analysis complete");
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Peptide,Source,LikelihoodRank,TPM,Transcripts");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writePeptideData(final PeptideData peptideData)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%.6f,%.6f,%s",
                    peptideData.Allele, peptideData.Peptide, peptideData.Source,
                    peptideData.LikelihoodRank, peptideData.mTpmTotal, peptideData.transcripts()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide output: {}", e.toString());
        }
    }

    private List<PeptideData> loadPeptideData(final String filename, final String source)
    {
        final List<PeptideData> peptideData = Lists.newArrayList();

        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("peptide search file({}) not found", filename);
            return peptideData;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            // Allele,Peptide,LikelihoodRank,GeneName,TransName - for proteome binders OR
            // Allele,Peptide,Source,LikelihoodRank - for binding data
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int rankIndex = fieldsIndexMap.get(FLD_LIKE_RANK);
            Integer transIndex = fieldsIndexMap.get("TransNames");

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIMITER, -1);

                String allele = values[alleleIndex];
                String peptide = values[peptideIndex];
                double likelihoodRank = Double.parseDouble(values[rankIndex]);

                List<String> transNames = Lists.newArrayList();

                if(transIndex != null)
                {
                    Arrays.stream(values[transIndex].split(ITEM_DELIM, -1)).forEach(x -> transNames.add(x));
                }

                peptideData.add(new PeptideData(allele, peptide, source, likelihoodRank, transNames));
            }

            NE_LOGGER.info("loaded {} peptides from file({})", peptideData.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read peptide data file: {}", e.toString());
        }

        return peptideData;
    }

    private Map<String,PeptideSearchData> loadPeptideSearchData(final String filename)
    {
        final Map<String,PeptideSearchData> peptideSearchDataMap = Maps.newHashMap();

        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("peptide search file({}) not found", filename);
            return peptideSearchDataMap;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            // Peptide,Genes,Transcripts,AminoAcidPos,UpFlank,DownFlank,MatchCount
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int transcriptsIndex = fieldsIndexMap.get("Transcripts");
            int genesIndex = fieldsIndexMap.get("Genes");

            for(String line :lines)
            {
                final String[] values = line.split(DELIMITER, -1);

                String peptide = values[peptideIndex];
                String transcriptsStr = values[transcriptsIndex];
                String genesStr = values[genesIndex];

                List<String> transcripts = Arrays.stream(transcriptsStr.split(ITEM_DELIM, -1)).collect(Collectors.toList());
                List<String> genes = Arrays.stream(genesStr.split(ITEM_DELIM, -1)).collect(Collectors.toList());

                peptideSearchDataMap.put(peptide, new PeptideSearchData(peptide, genes, transcripts));
            }

            NE_LOGGER.info("loaded {} peptide search records from file({})", peptideSearchDataMap.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read peptide search data file: {}", e.toString());
        }

        return peptideSearchDataMap;
    }

    private class PeptideData
    {
        public final String Allele;
        public final String Peptide;
        public final String Source;
        public final double LikelihoodRank;
        public final List<String> Transcripts;

        private double mTpmTotal;
        private boolean mHasTpm;

        public PeptideData(
                final String allele, final String peptide, final String source, final double likelihoodRank, final List<String> transcripts)
        {
            Allele = allele;
            Peptide = peptide;
            Source = source;
            LikelihoodRank = likelihoodRank;
            Transcripts = transcripts;
            mTpmTotal = 0;
            mHasTpm = false;
        }

        public void addTpm(final Double tpm)
        {
            if(tpm == null)
                return;

            mHasTpm = true;
            mTpmTotal += tpm;
        }

        public String toString() { return String.format("allele(%s) peptide(%s) source(%s)", Allele, Peptide, Source); }

        public String transcripts()
        {
            StringJoiner sj = new StringJoiner(ITEM_DELIM);
            Transcripts.forEach(x -> sj.add(x));
            return sj.toString();
        }
    }

    private class PeptideSearchData
    {
        public final String Peptide;
        public final List<String> Genes;
        public final List<String> Transcripts;

        public PeptideSearchData(final String peptide, final List<String> genes, final List<String> transcripts)
        {
            Peptide = peptide;
            Genes = genes;
            Transcripts = transcripts;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(IMMUNE_EXPRESSION_FILE, true, "Peptides to search for");
        options.addOption(PROTEOME_PEPTIDES_FILE, true, "Proteome binders file");
        options.addOption(VALIDATION_PEPTIDES_FILE, true, "Immunogenic peptide data file");
        options.addOption(PEPTIDE_SEARCH_FILE, true, "Peptide location in proteome");
        options.addOption(OUTPUT_ID, true, "Output file identifier");
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ExpressionAnalyser expressionAnalyser = new ExpressionAnalyser(cmd);
        expressionAnalyser.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
