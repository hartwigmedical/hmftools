package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.BIND_DELIM;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE_CFG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.neo.bind.TranscriptExpression;

import org.jetbrains.annotations.NotNull;

// routine for generating random peptides from the ref-genome for a given set of lengths and alleles (optional)
public class GenerateRandomPeptides
{
    private final int mRequiredPeptides;
    private final List<Integer> mPeptideLengths;
    private final List<String> mAlleles;

    private final EnsemblDataCache mEnsemblDataCache;
    private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;
    private final TranscriptExpression mTranscriptExpression;

    private final Set<String> mUniquePeptides;
    private final List<String> mPeptideExclusions;

    private final Random mRandom;
    private final int mMaxPeptidesPerGene;
    private final int mPeptideBaseGap;

    private int mPeptitdeLengthIndex;
    private boolean mAssignAllelesByFreq;

    private final String mOutputFile;
    private BufferedWriter mWriter;

    private static final String REQ_PEPTIDES = "req_peptides";
    private static final String PEPTIDES_LENGTHS = "peptide_lengths";
    private static final String ALLELES_FILE = "alleles_file";
    private static final String PEPTIDE_EXCLUSIONS_FILE = "peptide_exclusions_file";
    private static final String ASSIGN_ALLELES_BY_FREQ = "assign_alleles_by_freq";
    private static final String OUTPUT_FILE = "output_file";

    private static final int MAX_PER_GENE = 10;

    public GenerateRandomPeptides(final ConfigBuilder configBuilder)
    {
        String ensemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        mEnsemblDataCache = new EnsemblDataCache(ensemblDataDir, RefGenomeVersion.V37);
        mEnsemblDataCache.setRequiredData(false, false, false, true);

        mEnsemblDataCache.load(false);

        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList(), true);

        mTranscriptExpression = new TranscriptExpression(configBuilder.getValue(IMMUNE_EXPRESSION_FILE));

        mUniquePeptides = Sets.newHashSet();
        mPeptideExclusions = Lists.newArrayList();

        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mAssignAllelesByFreq = configBuilder.hasFlag(ASSIGN_ALLELES_BY_FREQ);

        mPeptitdeLengthIndex = 0;

        mRequiredPeptides = configBuilder.getInteger(REQ_PEPTIDES);

        mPeptideLengths = Lists.newArrayList();

        if(configBuilder.hasValue(PEPTIDES_LENGTHS))
        {
            Arrays.stream(configBuilder.getValue(PEPTIDES_LENGTHS).split(ITEM_DELIM, -1))
                    .forEach(x -> mPeptideLengths.add(Integer.parseInt(x)));
        }
        else
        {
            mPeptideLengths.add(9);
        }

        mAlleles = Lists.newArrayList();
        loadAlleleFrequencies(configBuilder.getValue(ALLELES_FILE));
        loadPeptideExclusions(configBuilder.getValue(PEPTIDE_EXCLUSIONS_FILE));

        mRandom = new Random();
        mPeptideBaseGap = MAX_PER_GENE - (int)round(log10(mRequiredPeptides)); // smaller gaps when more peptides are required
        mMaxPeptidesPerGene = (int)max(mRequiredPeptides / 20000.0 * 3, MAX_PER_GENE); // from about 20K coding canonical trans

        mWriter = initialiseWriter();
    }

    private void loadAlleleFrequencies(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(fileContents.get(0), BIND_DELIM);
            fileContents.remove(0);
            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            Integer freqIndex = fieldsIndexMap.get("AlleleFreq");

            for(String data : fileContents)
            {
                String[] items = data.split(BIND_DELIM);
                String allele = items[alleleIndex];

                mAlleles.add(allele);

                if(mAssignAllelesByFreq && freqIndex != null)
                {
                    int freq = Integer.parseInt(items[freqIndex]);

                    for(int i = 1; i < freq; ++i)
                    {
                        mAlleles.add(allele);
                    }
                }
            }

            NE_LOGGER.info("loaded {} alleles", fileContents.size());
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load alleles file: {}", e.toString());
        }
    }

    private void loadPeptideExclusions(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            mPeptideExclusions.addAll(Files.readAllLines(new File(filename).toPath()));
            mPeptideExclusions.remove(0);

            NE_LOGGER.info("loaded {} peptide exclusions", mPeptideExclusions.size());
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

        int reqPeptideLength = mPeptideLengths.stream().mapToInt(x -> x.intValue()).max().orElse(0);

        for(Map.Entry<String,List<GeneData>> entry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
        {
            String chromosome = entry.getKey();

            for(GeneData geneData : entry.getValue())
            {
                int genePeptideCount = 0;

                TranscriptData transData = mEnsemblDataCache.getTranscriptData(geneData.GeneId, "");

                if(transData == null || transData.CodingStart == null)
                    continue;

                if(mTranscriptExpression.hasData() && mTranscriptExpression.getExpression(transData.TransName) == null)
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
                int endPos = startPos + reqPeptideLength;

                while(endPos < codingLength)
                {
                    String peptide = aminoAcids.substring(startPos, endPos);

                    if(isValidPeptide(peptide))
                    {
                        mUniquePeptides.add(peptide);
                        writeData(chromosome, geneData.GeneName, startPos, endPos, peptide);

                        ++genePeptideCount;
                        ++totalPeptideCount;

                        if(totalPeptideCount > 0 && (totalPeptideCount % 10000) == 0)
                        {
                            NE_LOGGER.info("found {} peptides", totalPeptideCount);
                        }

                        if(genePeptideCount >= mMaxPeptidesPerGene || totalPeptideCount >= mRequiredPeptides)
                            break;
                    }

                    startPos = mPeptideBaseGap > 0 ? endPos + mRandom.nextInt(mPeptideBaseGap) : startPos + 1;
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

    private boolean isValidPeptide(final String peptide)
    {
        final String stopAminoAcid = String.valueOf(Codons.STOP_AMINO_ACID);

        if(peptide.contains(stopAminoAcid))
            return false;

        if(mUniquePeptides.contains(peptide))
            return false;

        for(int peptideLength : mPeptideLengths)
        {
            String peptideByLength = peptide.length() > peptideLength ? peptide.substring(0, peptideLength) : peptide;
            if(mPeptideExclusions.contains(peptideByLength))
                return false;
        }

        return true;
    }

    private String getAllele()
    {
        if(mAlleles.isEmpty())
            return "NONE";

        int index = mRandom.nextInt(mAlleles.size());
        return mAlleles.get(index);
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            if(!mAlleles.isEmpty())
                writer.write("Allele,");

            writer.write("Peptide,Chromosome,GeneName,PosStart,PosEnd");
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
            if(mAssignAllelesByFreq)
            {
                mWriter.write(String.format("%s,%s,%s,%s,%d,%d", getAllele(), peptide, chromosome, geneName, posStart, posEnd));
                mWriter.newLine();
            }
            else
            {
                if(mAlleles.isEmpty())
                {
                    mWriter.write(String.format("%s,%s,%s,%d,%d", peptide, chromosome, geneName, posStart, posEnd));
                    mWriter.newLine();
                }
                else
                {
                    for(String allele : mAlleles)
                    {
                        for(int peptideLength : mPeptideLengths)
                        {
                            String peptideByLength = peptide.length() > peptideLength ? peptide.substring(0, peptideLength) : peptide;

                            mWriter.write(String.format("%s,%s,%s,%s,%d,%d", allele, peptideByLength, chromosome, geneName, posStart, posEnd));
                            mWriter.newLine();
                        }
                    }
                }
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addEnsemblDir(configBuilder);
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addRequiredInteger(REQ_PEPTIDES, "Number of peptides to find randomly from the proteome");
        configBuilder.addConfigItem(PEPTIDES_LENGTHS, false, "Peptide lengths, separated by ';'");
        configBuilder.addPath(ALLELES_FILE, true, "File with alleles to assign");
        configBuilder.addPath(PEPTIDE_EXCLUSIONS_FILE, false, "File with training set peptides to avoid replicating");
        configBuilder.addPath(IMMUNE_EXPRESSION_FILE, false, IMMUNE_EXPRESSION_FILE_CFG);

        addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateRandomPeptides neoBinder = new GenerateRandomPeptides(configBuilder);
        neoBinder.run();
    }
}
