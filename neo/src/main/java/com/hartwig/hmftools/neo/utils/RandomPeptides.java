package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

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
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.neo.bind.BinderConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class RandomPeptides
{
    private final int mRequiredPeptides;
    private final List<Integer> mPeptideLengths;
    private final List<String> mAlleles;

    private final EnsemblDataCache mEnsemblDataCache;
    private final Map<String, TranscriptAminoAcids > mTransAminoAcidMap;

    private final Random mRandom;
    private final int mPeptideBaseGap;

    private int mPeptitdeLengthIndex;
    private int mAlleleIndex;

    private BufferedWriter mWriter;

    private static final String REQ_PEPTIDES = "req_peptides";
    private static final String PEPTIDES_LENGTHS = "peptide_lengths";
    private static final String ALLELES_FILE = "alleles_file";

    private static final int MAX_PER_GENE = 10;

    public RandomPeptides(final CommandLine cmd)
    {
        String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);
        mEnsemblDataCache = new EnsemblDataCache(ensemblDataDir, RefGenomeVersion.V37);
        mEnsemblDataCache.setRequiredData(true, false, false, true);

        mEnsemblDataCache.load(false);

        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList());

        mRandom = new Random();
        mPeptideBaseGap = 10;

        mPeptitdeLengthIndex = 0;
        mAlleleIndex = 0;

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

        try
        {
            String allelesFile = cmd.getOptionValue(ALLELES_FILE);
            final List<String> fileContents = Files.readAllLines(new File(allelesFile).toPath());
            fileContents.stream().filter(x -> !x.contains("Allele")).forEach(x -> mAlleles.add(x));

            NE_LOGGER.info("loaded {} alleles", mAlleles.size());
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load alleles file: {}", e.toString());
        }

        mWriter = initialiseWriter(parseOutputDir(cmd));
    }

    public void run()
    {
        NE_LOGGER.info("searching for {} peptides of lengths: {}", mRequiredPeptides, mPeptideLengths);

        int totalPeptideCount = 0;
        int transCodingCount = 0;
        int maxPerGene = (int)max(mRequiredPeptides / 25000.0 * 5, MAX_PER_GENE);

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

                    writeData(chromosome, geneData.GeneName, startPos, endPos, peptide);

                    ++genePeptideCount;
                    ++totalPeptideCount;

                    if(genePeptideCount >= maxPerGene || totalPeptideCount >= mRequiredPeptides)
                        break;

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
        if(mAlleleIndex >= mAlleles.size())
            mAlleleIndex = 0;

        String allele = mAlleles.get(mAlleleIndex);
        ++mAlleleIndex;
        return allele;
    }

    private static BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String filename = outputDir + "ref_genome_peptides.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

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
        options.addOption(REQ_PEPTIDES, true, "Output directory");
        options.addOption(PEPTIDES_LENGTHS, true, "Output directory");
        options.addOption(ALLELES_FILE, true, "File with alleles to assign");

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");


        BinderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        RandomPeptides neoBinder = new RandomPeptides(cmd);
        neoBinder.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
