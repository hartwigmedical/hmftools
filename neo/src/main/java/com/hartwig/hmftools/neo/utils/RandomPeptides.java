package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.epitope.NeoConfig.REF_GENOME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
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

    private final EnsemblDataCache mEnsemblDataCache;
    private final RefGenomeInterface mRefGenome;
    private final Random mRandom;

    private int mPeptitdeLengthIndex;

    private BufferedWriter mWriter;

    private static final String REQ_PEPTIDES = "req_peptides";
    private static final String PEPTIDES_LENGTHS = "peptide_lengths";
    private static final int MAX_PER_GENE = 10;

    public RandomPeptides(final CommandLine cmd)
    {
        mEnsemblDataCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), RefGenomeVersion.V37);
        mEnsemblDataCache.setRequiredData(true, false, false, true);

        mRefGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mRandom = new Random();
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

        mWriter = initialiseWriter(parseOutputDir(cmd));
    }

    private int getPeptideLength()
    {
        if(mPeptitdeLengthIndex >= mPeptideLengths.size())
            mPeptitdeLengthIndex = 0;

        int length = mPeptideLengths.get(mPeptitdeLengthIndex);
        ++mPeptitdeLengthIndex;
        return length;
    }

    public void run()
    {
        mEnsemblDataCache.load(false);

        NE_LOGGER.info("searching for {} peptides of lengths: {}", mRequiredPeptides, mPeptideLengths);

        int totalPeptideCount = 0;

        for(Map.Entry<String,List<GeneData>> entry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
        {
            String chromosome = entry.getKey();

            for(GeneData geneData : entry.getValue())
            {
                int genePeptideCount = 0;

                TranscriptData transData = mEnsemblDataCache.getTranscriptData(geneData.GeneId, "");

                if(transData == null || transData.CodingStart == null)
                    continue;

                for(ExonData exon : transData.exons())
                {
                    if(exon.Start < transData.CodingStart)
                        continue;

                    if(exon.End > transData.CodingEnd)
                        break;

                    int reqPeptideLength = getPeptideLength();
                    int requiredBases = reqPeptideLength * 3;

                    int maxEndPos = exon.End - requiredBases + 1;
                    int codingBases = maxEndPos - exon.Start;
                    if(codingBases <= 0)
                        continue;

                    int startPos = exon.Start + mRandom.nextInt(codingBases);
                    int endPos = startPos + requiredBases - 1;

                    String refCodingBases = mRefGenome.getBaseString(chromosome, startPos, endPos);

                    if(transData.posStrand())
                        refCodingBases = Nucleotides.reverseStrandBases(refCodingBases);

                    String aminoAcids = Codons.aminoAcidFromBases(refCodingBases);

                    writeData(chromosome, geneData.GeneName, startPos, endPos, aminoAcids);

                    ++genePeptideCount;
                    ++totalPeptideCount;

                    if(genePeptideCount >= MAX_PER_GENE || totalPeptideCount >= mRequiredPeptides)
                        break;
                }

                if(totalPeptideCount >= mRequiredPeptides)
                    break;
            }

            if(totalPeptideCount >= mRequiredPeptides)
                break;
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("found {} random peptides", totalPeptideCount);
    }

    private static BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String filename = outputDir + "ref_genome_peptides.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,GeneName,PosStart,PosEnd,Peptide");
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
            mWriter.write(String.format("%s,%s,%d,%d,%s", chromosome, geneName, posStart, posEnd, peptide));
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
        options.addOption(REF_GENOME, true, "Ref genome path");
        options.addOption(REQ_PEPTIDES, true, "Output directory");
        options.addOption(PEPTIDES_LENGTHS, true, "Output directory");
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
