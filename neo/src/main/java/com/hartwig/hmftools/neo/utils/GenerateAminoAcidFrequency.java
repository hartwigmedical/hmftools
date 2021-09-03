package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.neo.bind.BindConstants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

// class for generating amino acid frequencies from the proteome
public class GenerateAminoAcidFrequency
{
    private final Map<String,TranscriptAminoAcids > mTransAminoAcidMap;

    private final String mAminoAcidFreqFile;

    private final Map<Character,Double> mAminoAcidFrequencies;

    public static final String AMINO_ACID_FREQ_FILE = "amino_acid_freq_file";

    public GenerateAminoAcidFrequency(final CommandLine cmd)
    {
        mAminoAcidFrequencies = Maps.newHashMap();
        mAminoAcidFreqFile = cmd.getOptionValue(AMINO_ACID_FREQ_FILE);

        if(cmd.hasOption(ENSEMBL_DATA_DIR))
        {
            mTransAminoAcidMap = Maps.newHashMap();
            String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);
            EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList(), true);
        }
        else
        {
            mTransAminoAcidMap = null;
        }
    }

    public void generateFrequencies()
    {
        NE_LOGGER.info("measuring amino acid frequencies from full proteome: transcripts = {}", mTransAminoAcidMap.size());

        long totalAAs = 0;
        for(TranscriptAminoAcids transData : mTransAminoAcidMap.values())
        {
            for(int i = 0; i < transData.AminoAcids.length() - 1; ++i)
            {
                char aminoAcid = transData.AminoAcids.charAt(i);
                ++totalAAs;

                Double count = mAminoAcidFrequencies.get(aminoAcid);
                mAminoAcidFrequencies.put(aminoAcid, count != null ? count + 1 : 1);
            }
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(mAminoAcidFreqFile, false);

            writer.write("AminoAcid,Frequency,Percent");
            writer.newLine();

            for(Map.Entry<Character,Double> entry : mAminoAcidFrequencies.entrySet())
            {
                double freq = entry.getValue();
                writer.write(String.format("%c,%.0f,%.4f", entry.getKey(), freq, freq / (double)totalAAs));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise output file({}): {}", mAminoAcidFreqFile, e.toString());
        }

        NE_LOGGER.info("wrote amino acid frequencies");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addEnsemblDir(options);
        options.addOption(AMINO_ACID_FREQ_FILE, true, "Output filename");
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GenerateAminoAcidFrequency neoBinder = new GenerateAminoAcidFrequency(cmd);
        neoBinder.generateFrequencies();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
