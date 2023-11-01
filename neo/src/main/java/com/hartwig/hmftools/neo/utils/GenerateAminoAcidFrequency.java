package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

// class for generating amino acid frequencies from the proteome
public class GenerateAminoAcidFrequency
{
    private final Map<String,TranscriptAminoAcids > mTransAminoAcidMap;

    private final String mAminoAcidFreqFile;

    private final Map<Character,Double> mAminoAcidFrequencies;

    public static final String AMINO_ACID_FREQ_FILE = "amino_acid_freq_file";

    public GenerateAminoAcidFrequency(final ConfigBuilder configBuilder)
    {
        mAminoAcidFrequencies = Maps.newHashMap();
        mAminoAcidFreqFile = configBuilder.getValue(AMINO_ACID_FREQ_FILE);

        if(configBuilder.hasValue(ENSEMBL_DATA_DIR))
        {
            mTransAminoAcidMap = Maps.newHashMap();
            String ensemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addEnsemblDir(configBuilder);
        configBuilder.addPath(AMINO_ACID_FREQ_FILE, true, "Output filename");
        addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateAminoAcidFrequency neoBinder = new GenerateAminoAcidFrequency(configBuilder);
        neoBinder.generateFrequencies();
    }
}
