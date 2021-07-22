package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class AminoAcidFrequency
{
    private final Map<String,TranscriptAminoAcids > mTransAminoAcidMap;

    private final String mOutputFile;

    private final Map<Character,Integer> mAminoAcidFrequencies;

    private static final String OUTPUT_FILE = "output_file";

    public AminoAcidFrequency(final CommandLine cmd)
    {
        String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);

        mTransAminoAcidMap = Maps.newHashMap();
        mAminoAcidFrequencies = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList());

        mOutputFile = cmd.getOptionValue(OUTPUT_FILE);
    }

    public void run()
    {
        NE_LOGGER.info("measuring amino acid frequencies from full proteome: transcripts = {}", mTransAminoAcidMap.size());

        long totalAAs = 0;
        for(TranscriptAminoAcids transData : mTransAminoAcidMap.values())
        {
            for(int i = 0; i < transData.AminoAcids.length() - 1; ++i)
            {
                char aminoAcid = transData.AminoAcids.charAt(i);
                ++totalAAs;

                Integer count = mAminoAcidFrequencies.get(aminoAcid);
                mAminoAcidFrequencies.put(aminoAcid, count != null ? count + 1 : 1);
            }
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            writer.write("AminoAcid,Frequency,Percent");
            writer.newLine();

            for(Map.Entry<Character,Integer> entry : mAminoAcidFrequencies.entrySet())
            {
                int freq = entry.getValue();
                writer.write(String.format("%c,%d,%.4f", entry.getKey(), freq, freq / (double)totalAAs));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise output file({}): {}", mOutputFile, e.toString());
        }

        NE_LOGGER.info("wrote amino acid frequencies");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(ENSEMBL_DATA_DIR, true, "Ensembl data dir");
        options.addOption(OUTPUT_FILE, true, "Output filename");
        options.addOption(OUTPUT_DIR, true, "Output directory");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        AminoAcidFrequency neoBinder = new AminoAcidFrequency(cmd);
        neoBinder.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
