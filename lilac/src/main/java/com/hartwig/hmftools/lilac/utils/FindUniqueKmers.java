package com.hartwig.hmftools.lilac.utils;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_DEBUG;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_Y;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_REF_FILE;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.internal.$processor$.encode.$Eq;
import org.jetbrains.annotations.NotNull;

public class FindUniqueKmers
{
    private final String mResourceDir;

    private final List<HlaSequence> mAminoAcidSequences;
    private final List<HlaSequence> mHlaYAminoAcidSequences;

    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mHlaYNucleotideSequences;

    private final HlaAlleleCache mAlleleCache;

    private static final int KMER_MIN = 8;
    private static final int KMER_MAX = 20;

    public FindUniqueKmers(final CommandLine cmd)
    {
        mAlleleCache = new HlaAlleleCache();
        mResourceDir = cmd.getOptionValue(RESOURCE_DIR);

        mAminoAcidSequences = Lists.newArrayList();
        mHlaYAminoAcidSequences = Lists.newArrayList();

        mNucleotideSequences = Lists.newArrayList();
        mHlaYNucleotideSequences = Lists.newArrayList();
    }

    public void run()
    {
        // loadAminoAcidSequences();
        loadNucleotideSequences();

        LL_LOGGER.info("analysing {} HLA-Y alleles");

        if(!mAminoAcidSequences.isEmpty() && !mHlaYAminoAcidSequences.isEmpty())
        {
            LL_LOGGER.info("searching for unique k-mers");
            mHlaYAminoAcidSequences.forEach(x -> findUniqueKmers(x));
            LL_LOGGER.info("unique k-mer search complete");
        }

        if(!mNucleotideSequences.isEmpty() && !mHlaYNucleotideSequences.isEmpty())
        {
            LL_LOGGER.info("searching for unique loci");
            mHlaYNucleotideSequences.forEach(x -> findUniqueLoci(x));
            LL_LOGGER.info("unique loci search complete");
        }

        LL_LOGGER.info("HLA-Y analysis complete");
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

                if(!GENE_IDS.contains(allele.Gene) && !allele.Gene.equals(GENE_Y))
                    continue;

                String sequenceStr = items[1];

                HlaSequenceLoci newSequence = HlaSequenceFile.createFromReference(allele, sequenceStr, true);
                HlaSequence sequence = new HlaSequence(newSequence.Allele, newSequence.sequence());

                if(allele.Gene.equals(GENE_Y))
                    mHlaYAminoAcidSequences.add(sequence);
                else
                    mAminoAcidSequences.add(sequence);
            }
        }
        catch (IOException e)
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

                if(!GENE_IDS.contains(allele.Gene) && !allele.Gene.equals(GENE_Y))
                    continue;

                String sequenceStr = items[1];

                HlaSequenceLoci sequence = HlaSequenceFile.createFromReference(allele, sequenceStr, false);

                if(allele.Gene.equals(GENE_Y))
                    mHlaYNucleotideSequences.add(sequence);
                else
                    mNucleotideSequences.add(sequence);
            }
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", aminoAcidFilename, e.toString());
            return;
        }

        LL_LOGGER.info("loaded {} sequences from file {}", mAminoAcidSequences.size(), aminoAcidFilename);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("finding unique K-mer sequences");

        Options options = new Options();
        options.addOption(RESOURCE_DIR, true, "Path to resource files");
        options.addOption(OUTPUT_DIR, true, "Path to output");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        FindUniqueKmers findUniqueKmers = new FindUniqueKmers(cmd);
        findUniqueKmers.run();;


        LL_LOGGER.info("reference data written");
    }

}
