package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.OUTPUT_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConstants.ALL_NUCLEOTIDE_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.ALL_PROTEIN_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LOCI_POSITION;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
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
    private final LilacConfig mConfig;

    public final List<HlaSequenceLoci> NucleotideSequences;
    public final List<HlaSequenceLoci> AminoAcidSequences;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithInserts;
    public final List<HlaSequenceLoci> AminoAcidSequencesWithDeletes;

    private static final char SEQUENCE_DELIM = '|';
    private static final String NUC_REF_FILE = "lilac_ref_nucleotide_sequences.csv";
    private static final String AA_REF_FILE = "lilac_ref_aminoacid_sequences.csv";

    public ReferenceData(final String resourceDir, final LilacConfig config)
    {
        mResourceDir = resourceDir;
        mConfig = config;

        NucleotideSequences = Lists.newArrayList();
        AminoAcidSequences = Lists.newArrayList();
        AminoAcidSequencesWithInserts = Lists.newArrayList();
        AminoAcidSequencesWithDeletes = Lists.newArrayList();

        // load gene definitions and other constants
        populateHlaTranscripts();
        populateConstants();
    }

    public static void populateHlaTranscripts()
    {
        if(!HLA_TRANSCRIPTS.isEmpty())
            return;

        final Map<String, HmfTranscriptRegion> allTranscripts = HmfGenePanelSupplier.allGenesMap37();
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_A));
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_B));
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_C));
        LOCI_POSITION.initialise(HLA_TRANSCRIPTS);
    }

    public static void populateConstants()
    {
        for(Integer boundary : ALL_PROTEIN_EXON_BOUNDARIES)
        {
            ALL_NUCLEOTIDE_EXON_BOUNDARIES.add(boundary);
            ALL_NUCLEOTIDE_EXON_BOUNDARIES.add(boundary + 1);
            ALL_NUCLEOTIDE_EXON_BOUNDARIES.add(boundary + 2);
        }
    }

    public boolean load(boolean canUseConsolidated)
    {
        LL_LOGGER.info("Reading nucleotide files");

        String nucleotideFilename = mResourceDir + NUC_REF_FILE;

        if(canUseConsolidated && Files.exists(Paths.get(nucleotideFilename)))
        {
            if(!loadSequenceFile(nucleotideFilename, NucleotideSequences))
                return false;
        }
        else
        {
            NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/A_nuc.txt"));
            NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/B_nuc.txt"));
            NucleotideSequences.addAll(nucleotideLoci(mResourceDir + "/C_nuc.txt"));
        }

        LL_LOGGER.info("Reading protein files");

        String aminoAcidFilename = mResourceDir + AA_REF_FILE;

        if(canUseConsolidated && Files.exists(Paths.get(aminoAcidFilename)))
        {
            if(!loadSequenceFile(aminoAcidFilename, AminoAcidSequences))
                return false;
        }
        else
        {
            AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/A_prot.txt"));
            AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/B_prot.txt"));
            AminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "/C_prot.txt"));
        }

        AminoAcidSequencesWithInserts.addAll(AminoAcidSequences.stream().filter(x -> x.containsInserts()).collect(Collectors.toList()));
        AminoAcidSequencesWithDeletes.addAll(AminoAcidSequences.stream().filter(x -> x.containsDeletes()).collect(Collectors.toList()));

        return true;
    }

    private boolean excludeAllele(final HlaAllele allele)
    {
        final HlaAllele allele4d = allele.asFourDigit();

        if(EXCLUDED_ALLELES.stream().anyMatch(x -> allele4d.matches(x)))
            return true;

        if(mConfig == null)
            return false;

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            if(mConfig.ExpectedAlleles.stream().anyMatch(x -> x.matches(allele4d)))
                return false;

            if(mConfig.RestrictedAlleles.stream().noneMatch(x -> x.matches(allele4d)))
                return true;
        }

        return false;
    }

    private boolean loadSequenceFile(final String filename, final List<HlaSequenceLoci> sequenceData)
    {
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            fileContents.remove(0);

            for(String line : fileContents)
            {
                String[] items = line.split(",");

                if(items.length != 2)
                    return false;

                HlaAllele allele = HlaAllele.fromString(items[0]);

                if(excludeAllele(allele))
                    continue;

                List<String> sequences = Lists.newArrayList();
                String sequenceStr = items[1];

                int index = 0;
                String sequence = "";
                boolean inMulti = false;
                while(index < sequenceStr.length())
                {
                    char nextChar = sequenceStr.charAt(index);
                    boolean isMulti = nextChar == SEQUENCE_DELIM;

                    if(inMulti || isMulti)
                    {
                        if(inMulti && isMulti)
                        {
                            inMulti = false;
                            sequences.add(sequence);
                            sequence = "";
                        }
                        else if(isMulti)
                        {
                            // start of new multi-char sequence
                            inMulti = true;
                        }
                        else
                        {
                            sequence += nextChar;
                        }
                    }
                    else
                    {
                        sequences.add(String.valueOf(nextChar));
                    }

                    ++index;
                }

                sequenceData.add(new HlaSequenceLoci(allele, sequences));
            }

            LL_LOGGER.info("loaded {} from file {}", sequenceData.size(), filename);
            return true;
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to load ref sequence data from file({}): {}", filename, e.toString());
            return false;
        }
    }

    public List<HlaSequenceLoci> nucleotideLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .filter(x -> !excludeAllele(x.Allele))
                .collect(Collectors.toList());

        List<HlaSequence> reducedSequences = reduceToSixDigit(filteredSequences);

        return HlaSequenceLoci.create(reducedSequences);
    }

    public final List<HlaSequenceLoci> aminoAcidLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readFile(filename);

        List<HlaSequence> filteredSequences = sequences.stream()
                .filter(x -> !excludeAllele(x.Allele))
                //.filter(x -> EXCLUDED_ALLELES.stream().noneMatch(y -> x.Allele.asFourDigit().matches(y)))
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
        options.addOption(RESOURCE_DIR, true, "Path to resource files");
        options.addOption(OUTPUT_DIR, true, "Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        ReferenceData refData = new ReferenceData(cmd.getOptionValue(RESOURCE_DIR), null);

        EXCLUDED_ALLELES.clear();

        if(!refData.load(false))
        {
            LL_LOGGER.error("ref data loading failed");
            System.exit(1);
        }

        final String outputDir = parseOutputDir(cmd);
        refData.rewriteRefData(outputDir);

        LL_LOGGER.info("Reference data written");
    }

    public void rewriteRefData(final String outputDir)
    {
        String nucleotideSequenceFile = outputDir + "lilac_ref_nucleotide_sequences.csv";
        writeSequenceData(nucleotideSequenceFile, NucleotideSequences);

        String aminoAcidSequenceFile = outputDir + "lilac_ref_aminoacid_sequences.csv";
        writeSequenceData(aminoAcidSequenceFile, AminoAcidSequences);
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
                writer.write(sequenceLoci.getAllele() + ",");

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

