package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.OUTPUT_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConstants.EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.SEQUENCE_DELIM;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit;

import java.io.BufferedWriter;
import java.io.IOException;
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

public class GenerateReferenceSequences
{
    private final String mResourceDir;

    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mAminoAcidSequences;

    private final HlaAlleleCache mAlleleCache;

    public GenerateReferenceSequences(final String resourceDir)
    {
        mAlleleCache = new HlaAlleleCache();
        mResourceDir = resourceDir;

        mNucleotideSequences = Lists.newArrayList();
        mAminoAcidSequences = Lists.newArrayList();

    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("generating HLA sequence reference data");

        Options options = new Options();
        options.addOption(RESOURCE_DIR, true, "Path to resource files");
        options.addOption(OUTPUT_DIR, true, "Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        GenerateReferenceSequences seqGenerator = new GenerateReferenceSequences(cmd.getOptionValue(RESOURCE_DIR));

        EXCLUDED_ALLELES.clear();

        if(!seqGenerator.loadSequenceFiles())
        {
            LL_LOGGER.error("ref data loading failed");
            System.exit(1);
        }

        final String outputDir = parseOutputDir(cmd);
        seqGenerator.rewriteRefData(outputDir);

        LL_LOGGER.info("reference data written");
    }

    public void rewriteRefData(final String outputDir)
    {
        String nucleotideSequenceFile = outputDir + NUC_REF_FILE;
        writeSequenceData(nucleotideSequenceFile, mNucleotideSequences);

        String aminoAcidSequenceFile = outputDir + AA_REF_FILE;
        writeSequenceData(aminoAcidSequenceFile, mAminoAcidSequences);
    }

    public boolean loadSequenceFiles()
    {
        LL_LOGGER.info("reading nucleotide files");

        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/A_nuc.txt"));
        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/B_nuc.txt"));
        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/C_nuc.txt"));
        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/Y_nuc.txt"));

        LL_LOGGER.info("reading protein files");

        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/A_prot.txt"));
        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/B_prot.txt"));
        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/C_prot.txt"));

        return true;
    }

    private List<HlaSequenceLoci> nucleotideLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .collect(Collectors.toList());

        List<HlaSequence> reducedSequences = reduceToSixDigit(filteredSequences);

        return buildSequences(reducedSequences, false);
    }

    private List<HlaSequenceLoci> aminoAcidLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .collect(Collectors.toList());

        List<HlaSequence> reducedSequences = reduceToFourDigit(filteredSequences);

        reducedSequences = reducedSequences.stream()
                .map(x -> x.getRawSequence().endsWith("X") ? x : x.copyWithAdditionalSequence("X"))
                .collect(Collectors.toList());

        return buildSequences(reducedSequences, true);
    }

    private List<HlaSequenceLoci> buildSequences(final List<HlaSequence> sequences, boolean isProteinFile)
    {
        // the first A, B and C entries in both the nucleotide and amino acid files (A*01:01:01, B*07:02:01, C*01:02:01)
        // are used as a reference for all the others
        final String reference = sequences.get(0).getRawSequence();

        // build sequences using the allele cache to avoid any duplication of alleles
        List<HlaSequenceLoci> newSequences = Lists.newArrayList();

        for(HlaSequence sequence : sequences)
        {
            HlaAllele allele = isProteinFile ?
                    mAlleleCache.requestFourDigit(sequence.Allele.toString()) : mAlleleCache.request(sequence.Allele.toString());

            HlaSequenceLoci newSequence = HlaSequenceLoci.create(allele, sequence.getRawSequence(), reference);

            if(!newSequence.getSequences().isEmpty())
            {
                newSequences.add(newSequence);
            }
        }

        return newSequences;
    }

    private void writeSequenceData(final String filename, final List<HlaSequenceLoci> sequenceData)
    {
        try
        {
            final BufferedWriter writer = createBufferedWriter(filename, false);

            // header
            writer.write("Allele,Sequence");
            writer.newLine();

            for(HlaSequenceLoci sequenceLoci : sequenceData)
            {
                writer.write(sequenceLoci.Allele + ",");

                for(String sequence : sequenceLoci.getSequences())
                {
                    if(sequence.length() == 1)
                    {
                        writer.write(sequence);
                    }
                    else
                    {
                        writer.write(SEQUENCE_DELIM + sequence + SEQUENCE_DELIM);
                    }
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write file {}: {}", filename, e.toString());
        }
    }

}
