package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.OUTPUT_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConstants.EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class ReferenceData
{
    private final String mResourceDir;

    public final List<HlaSequenceLoci> NucleotideSequences;
    public final List<HlaSequenceLoci> AminoAcidSequences;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithInserts;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithDeletes;

    public ReferenceData(final String resourceDir)
    {
        mResourceDir = resourceDir;

        NucleotideSequences = Lists.newArrayList();
        AminoAcidSequences = Lists.newArrayList();
        AminoAcidSequencesWithInserts = Lists.newArrayList();
        AminoAcidSequencesWithDeletes = Lists.newArrayList();
    }

    public boolean load()
    {
        LL_LOGGER.info("Reading nucleotide files");
        NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/A_nuc.txt"));
        NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/B_nuc.txt"));
        NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/C_nuc.txt"));

        LL_LOGGER.info("Reading protein files");
        AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/A_prot.txt"));
        AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/B_prot.txt"));
        AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/C_prot.txt"));

        AminoAcidSequencesWithInserts.addAll(AminoAcidSequences.stream().filter(x -> x.containsInserts()).collect(Collectors.toList()));
        AminoAcidSequencesWithDeletes.addAll(AminoAcidSequences.stream().filter(x -> x.containsDeletes()).collect(Collectors.toList()));

        return true;
    }

    public final List<HlaSequenceLoci> nucleotideLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .filter(x -> EXCLUDED_ALLELES.stream().noneMatch(y -> x.Allele.asFourDigit().matches(y)))
                .collect(Collectors.toList());

        List<HlaSequence> reducedSequences = reduceToSixDigit(filteredSequences);

        return HlaSequenceLoci.create(reducedSequences);
    }

    public static final List<HlaSequenceLoci> aminoAcidLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .filter(x -> EXCLUDED_ALLELES.stream().noneMatch(y -> x.Allele.asFourDigit().matches(y)))
                .collect(Collectors.toList());

        List<HlaSequence> reducedSequences = reduceToFourDigit(filteredSequences);

        reducedSequences = reducedSequences.stream()
                .map(x -> x.getRawSequence().endsWith("X") ? x : x.copyWithAdditionalSequence("X"))
                .collect(Collectors.toList());

        return HlaSequenceLoci.create(reducedSequences);
    }

    // reference data rewrite
    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("Rewriting reference data");

        Options options = new Options();
        options.addOption(RESOURCE_DIR, true,"Path to resource files");
        options.addOption(OUTPUT_DIR, true,"Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        ReferenceData refData = new ReferenceData(cmd.getOptionValue(RESOURCE_DIR));
        refData.rewriteRefData();

    }

    public void rewriteRefData()
    {
        // public final List<HlaSequenceLoci> NucleotideSequences;
        // public final List<HlaSequenceLoci> AminoAcidSequences;

    }

}

