package com.hartwig.hmftools.lilac.utils;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.MHC_CLASS;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.registerCommonConfig;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.loadHlaTranscripts;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_Y;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.lilac.GeneCache;
import com.hartwig.hmftools.lilac.MhcClass;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.jetbrains.annotations.NotNull;

public class FindUniqueKmers
{
    private final String mResourceDir;
    private final MhcClass mClassType;

    private final GeneCache mGeneCache;

    private final List<HlaSequence> mAminoAcidSequences;
    private final List<HlaSequence> mHlaYAminoAcidSequences;

    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mHlaYNucleotideSequences;

    private final HlaAlleleCache mAlleleCache;

    private static final int KMER_MIN = 8;
    private static final int KMER_MAX = 20;

    private FindUniqueKmers(final ConfigBuilder configBuilder)
    {
        mAlleleCache = new HlaAlleleCache();
        mResourceDir = configBuilder.getValue(RESOURCE_DIR);

        mClassType = MhcClass.valueOf(configBuilder.getValue(MHC_CLASS));

        Map<HlaGene, TranscriptData> hlaTranscriptMap = loadHlaTranscripts(RefGenomeVersion.V37, mClassType);
        mGeneCache = new GeneCache(mClassType, hlaTranscriptMap);

        mAminoAcidSequences = Lists.newArrayList();
        mHlaYAminoAcidSequences = Lists.newArrayList();

        mNucleotideSequences = Lists.newArrayList();
        mHlaYNucleotideSequences = Lists.newArrayList();
    }

    private void run()
    {
        LL_LOGGER.info("finding unique K-mer sequences");

        loadNucleotideSequences();

        LL_LOGGER.info("analysing {} HLA-Y alleles");

        if(!mAminoAcidSequences.isEmpty() && !mHlaYAminoAcidSequences.isEmpty())
        {
            LL_LOGGER.info("searching for unique k-mers");
            mHlaYAminoAcidSequences.forEach(this::findUniqueKmers);
            LL_LOGGER.info("unique k-mer search complete");
        }

        if(!mNucleotideSequences.isEmpty() && !mHlaYNucleotideSequences.isEmpty())
        {
            LL_LOGGER.info("searching for unique loci");
            mHlaYNucleotideSequences.forEach(this::findUniqueLoci);
            LL_LOGGER.info("unique loci search complete");
        }

        LL_LOGGER.info("search complete");
    }

    private void findUniqueKmers(final HlaSequence sequence)
    {
        LL_LOGGER.info("analysing allele({})", sequence.Allele);

        String sequenceStr = sequence.getRawSequence();
        int sequenceLength = sequenceStr.length();

        for(int length = KMER_MIN; length <= KMER_MAX; ++length)
        {
            for(int i = 0; i < sequenceLength; ++i)
            {
                if(i + length >= sequenceLength)
                    break;

                String kmer = sequenceStr.substring(i, i + length);
                boolean matched = searchKmer(kmer);

                if(!matched)
                {
                    LL_LOGGER.debug("HLA-Y allele({}) unmatched kmer({}) loci({} - {})",
                            sequence.Allele, kmer, i, i + length - 1);
                }
            }
        }
    }

    private boolean searchKmer(final String kmer)
    {
        boolean matched = false;

        for(HlaSequence sequence : mAminoAcidSequences)
        {
            if(sequence.getRawSequence().contains(kmer))
            {
                matched = true;
                break;
            }
        }

        return matched;
    }

    private void findUniqueLoci(final HlaSequenceLoci hlaySequence)
    {
        for(int locus = 0; locus < hlaySequence.length(); ++locus)
        {
            String nucleotide = hlaySequence.sequence(locus);

            if(nucleotide.equals(DEL_STR))
                continue;

            boolean inOtherHlaY = false;
            boolean inOtherAllele = false;

            for(HlaSequenceLoci sequence : mNucleotideSequences)
            {
                if(sequence.length() <= locus)
                    continue;

                if(sequence.sequence(locus).equals(nucleotide))
                {
                    inOtherAllele = true;
                    break;
                }
            }

            for(HlaSequenceLoci sequence : mHlaYNucleotideSequences)
            {
                if(sequence == hlaySequence)
                    continue;

                if(sequence.length() <= locus)
                    continue;

                if(sequence.sequence(locus).equals(nucleotide))
                {
                    inOtherHlaY = true;
                    break;
                }
            }

            if(!inOtherAllele && !inOtherHlaY)
            {
                LL_LOGGER.debug("HLA-Y allele({}) locus({}) unique vs other({}) hlaY({})",
                        hlaySequence.Allele, locus, inOtherAllele, inOtherHlaY);
            }
        }
    }

    private void loadAminoAcidSequences()
    {
        String aminoAcidFilename = mResourceDir + AA_REF_FILE;

        LL_LOGGER.info("loading protein file: {}", aminoAcidFilename);

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(aminoAcidFilename).toPath());
            fileContents.remove(0); // remove header

            for(String line : fileContents)
            {
                String[] items = line.split(",");

                if(items.length != 2)
                    return;

                String alleleStr = items[0];

                HlaAllele allele = mAlleleCache.requestFourDigit(alleleStr);

                if(!mGeneCache.GeneNames.contains(allele.Gene) && allele.Gene != HLA_Y)
                    continue;

                String sequenceStr = items[1];

                HlaSequenceLoci newSequence = HlaSequenceFile.createFromReference(allele, sequenceStr, true);
                HlaSequence sequence = new HlaSequence(newSequence.Allele, newSequence.sequence());

                if(allele.Gene == HLA_Y)
                    mHlaYAminoAcidSequences.add(sequence);
                else
                    mAminoAcidSequences.add(sequence);
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", aminoAcidFilename, e.toString());
            return;
        }

        LL_LOGGER.info("loaded {} sequences from file {}", mAminoAcidSequences.size(), aminoAcidFilename);
    }

    private void loadNucleotideSequences()
    {
        String aminoAcidFilename = mResourceDir + NUC_REF_FILE;

        LL_LOGGER.info("loading nucleotide file: {}", aminoAcidFilename);

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(aminoAcidFilename).toPath());
            fileContents.remove(0); // remove header

            for(String line : fileContents)
            {
                String[] items = line.split(",");

                if(items.length != 2)
                    return;

                String alleleStr = items[0];

                HlaAllele allele = mAlleleCache.request(alleleStr);

                if(!mGeneCache.GeneNames.contains(allele.Gene) && allele.Gene != HLA_Y)
                    continue;

                String sequenceStr = items[1];

                HlaSequenceLoci sequence = HlaSequenceFile.createFromReference(allele, sequenceStr, false);

                if(allele.Gene == HLA_Y)
                    mHlaYNucleotideSequences.add(sequence);
                else
                    mNucleotideSequences.add(sequence);
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", aminoAcidFilename, e.toString());
            return;
        }

        LL_LOGGER.info("loaded {} sequences from file {}", mAminoAcidSequences.size(), aminoAcidFilename);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        registerCommonConfig(configBuilder);

        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FindUniqueKmers findUniqueKmers = new FindUniqueKmers(configBuilder);
        findUniqueKmers.run();
    }

}
