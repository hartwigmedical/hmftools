package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_LIKE_RANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE_CFG;
import static com.hartwig.hmftools.neo.utils.PeptideExpressionData.SOURCE_PROTEOME;
import static com.hartwig.hmftools.neo.utils.PeptideExpressionData.SOURCE_VALIDATION;

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
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.BindCommon;
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
    private final List<PeptideExpressionData> mProteomePeptides;
    private final List<PeptideExpressionData> mBinderPeptides;
    private final ExpressionDistribution mExpressionDistribution;
    private final GenePropensity mGenePropensity;

    private final Set<String> mValidationAlleles;

    private final String mOutputDir;
    private final String mOutputId;

    private final BufferedWriter mWriter;

    private static final String PROTEOME_PEPTIDES_FILE = "proteome_binders_file";
    private static final String VALIDATION_PEPTIDES_FILE = "validation_file";
    private static final String PEPTIDE_SEARCH_FILE = "peptide_search_file";

    public ExpressionAnalyser(final CommandLine cmd)
    {
        mTranscriptExpression = new TranscriptExpression(cmd.getOptionValue(IMMUNE_EXPRESSION_FILE));

        mBinderPeptides = loadPeptideData(cmd.getOptionValue(VALIDATION_PEPTIDES_FILE), SOURCE_VALIDATION);
        mValidationAlleles = mBinderPeptides.stream().map(x -> x.Allele).collect(Collectors.toSet());

        mProteomePeptides = loadPeptideData(cmd.getOptionValue(PROTEOME_PEPTIDES_FILE), SOURCE_PROTEOME);

        mPeptideSearchDataMap = loadPeptideSearchData(cmd.getOptionValue(PEPTIDE_SEARCH_FILE));

        mExpressionDistribution = new ExpressionDistribution();

        mGenePropensity = new GenePropensity(mExpressionDistribution, cmd.getOptionValue(ENSEMBL_DATA_DIR));

        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        String outputFile = BindCommon.formFilename(mOutputDir, "peptide_expression", mOutputId);
        mWriter = initialiseWriter(outputFile);
    }

    public void run()
    {
        if(!mTranscriptExpression.hasData() || mPeptideSearchDataMap.isEmpty() || mProteomePeptides.isEmpty() || mBinderPeptides.isEmpty())
        {
            NE_LOGGER.error("failed to initialise");
            System.exit(1);
        }

        NE_LOGGER.info("calculating expression for {} binding peptides from {} alleles",
                mBinderPeptides.size(), mValidationAlleles.size());

        for(PeptideExpressionData pepExpData : mBinderPeptides)
        {
            // find transcripts based on search results
            PeptideSearchData searchData = mPeptideSearchDataMap.get(pepExpData.Peptide);
            if(searchData == null)
            {
                NE_LOGGER.warn("peptide({}) missing search data");
                continue;
            }

            pepExpData.Transcripts.addAll(searchData.Transcripts);

            for(String transName : pepExpData.Transcripts)
            {
                pepExpData.addTpm(transName, mTranscriptExpression.getExpression(transName));
            }

            mExpressionDistribution.process(pepExpData);
            writePeptideData(pepExpData);
        }

        NE_LOGGER.info("calculating expression for {} proteome peptides", mProteomePeptides.size());

        for(int i = 0; i < mProteomePeptides.size(); ++i)
        {
            PeptideExpressionData pepExpData = mProteomePeptides.get(i);

            if(!mValidationAlleles.contains(pepExpData.Allele))
                continue;

            for(String transName : pepExpData.Transcripts)
            {
                pepExpData.addTpm(transName, mTranscriptExpression.getExpression(transName));
            }

            mExpressionDistribution.process(pepExpData);
            writePeptideData(pepExpData);

            if(i > 0 && (i % 100000) == 0)
            {
                NE_LOGGER.info("processed {} proteome peptides", i);
            }
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("generating and writing expression distributions");

        mExpressionDistribution.formDistributions();
        mExpressionDistribution.writeDistributions(mOutputDir, mOutputId);

        NE_LOGGER.info("calculating gene propensity");

        mBinderPeptides.forEach(x -> mGenePropensity.process(x));
        mProteomePeptides.forEach(x -> mGenePropensity.process(x));

        mGenePropensity.writeResults(mOutputDir, mOutputId);

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

    private void writePeptideData(final PeptideExpressionData peptideData)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%.6f,%.6f,%s",
                    peptideData.Allele, peptideData.Peptide, peptideData.Source,
                    peptideData.LikelihoodRank, peptideData.hasTpm() ? peptideData.tpm() : -1, peptideData.transcripts()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide output: {}", e.toString());
        }
    }

    private List<PeptideExpressionData> loadPeptideData(final String filename, final String source)
    {
        final List<PeptideExpressionData> peptideData = Lists.newArrayList();

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

                peptideData.add(new PeptideExpressionData(allele, peptide, source, likelihoodRank, transNames));
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
        options.addOption(IMMUNE_EXPRESSION_FILE, true, IMMUNE_EXPRESSION_FILE_CFG);
        options.addOption(PROTEOME_PEPTIDES_FILE, true, "Proteome binders file");
        options.addOption(VALIDATION_PEPTIDES_FILE, true, "Immunogenic peptide data file");
        options.addOption(PEPTIDE_SEARCH_FILE, true, "Peptide location in proteome");
        options.addOption(OUTPUT_ID, true, "Output file identifier");
        addEnsemblDir(options);
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
