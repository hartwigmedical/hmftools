package com.hartwig.hmftools.lilac.utils;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_Y;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;

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
import org.jetbrains.annotations.NotNull;

public class FindUniqueKmers
{
    private final String mResourceDir;

    private final List<HlaSequence> mAminoAcidSequences;
    private final List<HlaSequence> mHlaYSequences;

    private final HlaAlleleCache mAlleleCache;

    private static final int KMER_MIN = 8;
    private static final int KMER_MAX = 20;

    public FindUniqueKmers(final CommandLine cmd)
    {
        mAlleleCache = new HlaAlleleCache();
        mResourceDir = cmd.getOptionValue(RESOURCE_DIR);

        mAminoAcidSequences = Lists.newArrayList();
        mHlaYSequences = Lists.newArrayList();
    }

    public void run()
    {
        loadAminoAcidSequences();

        if(mAminoAcidSequences.isEmpty())
            return;

        mHlaYSequences.addAll(mAminoAcidSequences.stream().filter(x -> x.Allele.Gene.equals(GENE_Y)).collect(Collectors.toList()));
        mHlaYSequences.forEach(x -> mAminoAcidSequences.remove(x));

        LL_LOGGER.info("analysing {} HLA-Y alleles");

        mHlaYSequences.forEach(x -> analyseAllele(x));

        LL_LOGGER.info("HLA-Y kmer search complete");
    }

    private void analyseAllele(final HlaSequence sequence)
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
                    LL_LOGGER.info("HLA-Y allele({}) unmatched kmer({}) loci({} - {})",
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

                String sequenceStr = items[1];

                HlaSequenceLoci newSequence = HlaSequenceFile.createFromReference(allele, sequenceStr, true);
                HlaSequence sequence = new HlaSequence(newSequence.Allele, newSequence.sequence());
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("finding unique K-mer sequences");

        Options options = new Options();
        options.addOption(RESOURCE_DIR, true, "Path to resource files");
        options.addOption(OUTPUT_DIR, true, "Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        FindUniqueKmers findUniqueKmers = new FindUniqueKmers(cmd);
        findUniqueKmers.run();;


        LL_LOGGER.info("reference data written");
    }

}
