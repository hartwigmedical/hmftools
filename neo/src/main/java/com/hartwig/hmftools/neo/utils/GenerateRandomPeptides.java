package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class GenerateRandomPeptides
{
    private final int mRequiredPeptides;
    private final List<Integer> mPeptideLengths;
    private final List<String> mAlleles;

    private final EnsemblDataCache mEnsemblDataCache;
    private final Map<String,TranscriptAminoAcids > mTransAminoAcidMap;

    private final Random mRandom;
    private final int mMaxPeptidesPerGene;
    private final int mPeptideBaseGap;

    private int mPeptitdeLengthIndex;
    private int mAlleleIndex;

    private final String mOutputFile;
    private BufferedWriter mWriter;

    private static final String REQ_PEPTIDES = "req_peptides";
    private static final String PEPTIDES_LENGTHS = "peptide_lengths";
    private static final String ALLELES_FILE = "alleles_file";
    private static final String OUTPUT_FILE = "output_file";

    private static final int MAX_PER_GENE = 10;

    public GenerateRandomPeptides(final CommandLine cmd)
    {
        String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);
        mEnsemblDataCache = new EnsemblDataCache(ensemblDataDir, RefGenomeVersion.V37);
        mEnsemblDataCache.setRequiredData(true, false, false, true);

        mEnsemblDataCache.load(false);

        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList());

        mOutputFile = cmd.getOptionValue(OUTPUT_FILE);

        mPeptitdeLengthIndex = 0;

        mRequiredPeptides = Integer.parseInt(cmd.getOptionValue(REQ_PEPTIDES));

        mPeptideLengths = Lists.newArrayList();

        if(cmd.hasOption(PEPTIDES_LENGTHS))
        {
            Arrays.stream(cmd.getOptionValue(PEPTIDES_LENGTHS).split(";", -1))
                    .forEach(x -> mPeptideLengths.add(Integer.parseInt(x)));
        }
        else
        {
            mPeptideLengths.add(9);
        }

        mAlleles = Lists.newArrayList();
        loadAlleleFrequencies(cmd.getOptionValue(ALLELES_FILE));

        mRandom = new Random();
        mPeptideBaseGap = MAX_PER_GENE - (int)round(log10(mRequiredPeptides)); // smaller gaps when more peptides are required
        mMaxPeptidesPerGene = (int)max(mRequiredPeptides / 20000.0 * 3, MAX_PER_GENE); // from about 20K coding canonical trans

        mWriter = initialiseWriter(parseOutputDir(cmd));
    }

    private void loadAlleleFrequencies(final String filename)
    {
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(fileContents.get(0), DELIM);
            fileContents.remove(0);
            int alleleIndex = fieldsIndexMap.get("Allele");
            int freqIndex = fieldsIndexMap.get("AlleleFreq");

            for(String data : fileContents)
            {
                String[] items = data.split(DELIM);
                String allele = items[alleleIndex];
                int freq = Integer.parseInt(items[freqIndex]);

                for(int i = 0; i < freq; ++i)
                {
                    mAlleles.add(allele);
                }
            }

            NE_LOGGER.info("loaded {} alleles", fileContents.size());
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load alleles file: {}", e.toString());
        }
    }

    public void run()
    {
        NE_LOGGER.info("searching for {} peptides of lengths: {}", mRequiredPeptides, mPeptideLengths);

        int totalPeptideCount = 0;
        int transCodingCount = 0;

        final String stopAminoAcid = String.valueOf(Codons.STOP_AMINO_ACID);

        for(Map.Entry<String,List<GeneData>> entry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
        {
            String chromosome = entry.getKey();

            for(GeneData geneData : entry.getValue())
            {
                int genePeptideCount = 0;

                TranscriptData transData = mEnsemblDataCache.getTranscriptData(geneData.GeneId, "");

                if(transData == null || transData.CodingStart == null)
                    continue;

                ++transCodingCount;

                TranscriptAminoAcids transAminoAcids = mTransAminoAcidMap.get(transData.TransName);

                if(transAminoAcids == null)
                {
                    NE_LOGGER.error("gene({}) trans({}) missing amino-acid sequence", geneData.GeneId, transData.TransName);
                    continue;
                }

                final String aminoAcids = transAminoAcids.AminoAcids;
                int codingLength = aminoAcids.length() - 1; // excluding the stop codon

                int startPos = mPeptideBaseGap > 0 ? mRandom.nextInt(mPeptideBaseGap) : 0;
                int reqPeptideLength = getPeptideLength();
                int endPos = startPos + reqPeptideLength;

                while(endPos < codingLength)
                {
                    String peptide = aminoAcids.substring(startPos, endPos);

                    if(!peptide.contains(stopAminoAcid))
                    {
                        writeData(chromosome, geneData.GeneName, startPos, endPos, peptide);

                        ++genePeptideCount;
                        ++totalPeptideCount;

                        if(genePeptideCount >= mMaxPeptidesPerGene || totalPeptideCount >= mRequiredPeptides)
                            break;
                    }

                    startPos = mPeptideBaseGap > 0 ? endPos + mRandom.nextInt(mPeptideBaseGap) : startPos + 1;
                    reqPeptideLength = getPeptideLength();
                    endPos = startPos + reqPeptideLength;
                }

                if(totalPeptideCount >= mRequiredPeptides)
                    break;
            }

            if(totalPeptideCount >= mRequiredPeptides)
                break;
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("found {} random peptides from {} coding transcripts", totalPeptideCount, transCodingCount);
    }

    private int getPeptideLength()
    {
        if(mPeptitdeLengthIndex >= mPeptideLengths.size())
            mPeptitdeLengthIndex = 0;

        int length = mPeptideLengths.get(mPeptitdeLengthIndex);
        ++mPeptitdeLengthIndex;
        return length;
    }

    private String getAllele()
    {
        int index = mRandom.nextInt(mAlleles.size());
        return mAlleles.get(index);

        /*
        if(mAlleleIndex >= mAlleles.size())
            mAlleleIndex = 0;

        String allele = mAlleles.get(mAlleleIndex);
        ++mAlleleIndex;
        return allele;
        */
    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            writer.write("Allele,Peptide,Chromosome,GeneName,PosStart,PosEnd");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writeData(final String chromosome, final String geneName, int posStart, int posEnd, final String peptide)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                    getAllele(), peptide, chromosome, geneName, posStart, posEnd));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(ENSEMBL_DATA_DIR, true, "Ensembl data dir");
        options.addOption(OUTPUT_FILE, true, "Output filename");
        options.addOption(REQ_PEPTIDES, true, "Number of peptides to find randomly from the proteome");
        options.addOption(PEPTIDES_LENGTHS, true, "Peptide lengths, separated by ';'");
        options.addOption(ALLELES_FILE, true, "File with alleles to assign");

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GenerateRandomPeptides neoBinder = new GenerateRandomPeptides(cmd);
        neoBinder.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
