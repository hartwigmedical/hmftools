package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.neo.bind.BindConstants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

// class for generating then reading back in amino acid frequencies from the proteome
public class AminoAcidFrequency
{
    private final Map<String,TranscriptAminoAcids > mTransAminoAcidMap;

    private final String mAminoAcidFreqFile;

    private final Map<Character,Double> mAminoAcidFrequencies;

    public static final String AMINO_ACID_FREQ_FILE = "amino_acid_freq_file";

    public AminoAcidFrequency(final CommandLine cmd)
    {
        mAminoAcidFrequencies = Maps.newHashMap();
        mAminoAcidFreqFile = cmd.getOptionValue(AMINO_ACID_FREQ_FILE);

        if(cmd.hasOption(ENSEMBL_DATA_DIR))
        {
            mTransAminoAcidMap = Maps.newHashMap();
            String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);
            EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList());
        }
        else
        {
            mTransAminoAcidMap = null;
        }
    }

    public Map<Character,Double> getAminoAcidFrequencies() { return mAminoAcidFrequencies; }

    public double getAminoAcidFrequency(final char aminoAcid)
    {
        Double percent = mAminoAcidFrequencies.get(aminoAcid);
        return percent != null ? percent : 1.0 / BindConstants.AMINO_ACIDS.size();
    }

    public void loadFrequencies()
    {
        if(mAminoAcidFreqFile == null)
            return;

        try
        {
            final List<String> lines = Files.readAllLines(new File(mAminoAcidFreqFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int aaIndex = fieldsIndexMap.get("AminoAcid");
            int percentIndex = fieldsIndexMap.get("Percent");

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                char aminoAcid = items[aaIndex].charAt(0);
                double percent = Double.parseDouble(items[percentIndex]);
                mAminoAcidFrequencies.put(aminoAcid, percent);
            }

            NE_LOGGER.info("loaded {} amino-acid frequencies from file({})",
                    mAminoAcidFrequencies.size(), mAminoAcidFreqFile);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load amino-acid frequencies from file: {}", e.toString());
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
        options.addOption(ENSEMBL_DATA_DIR, true, "Ensembl data dir");
        options.addOption(AMINO_ACID_FREQ_FILE, true, "Output filename");
        options.addOption(OUTPUT_DIR, true, "Output directory");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        AminoAcidFrequency neoBinder = new AminoAcidFrequency(cmd);
        neoBinder.generateFrequencies();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
