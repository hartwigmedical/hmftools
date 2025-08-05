package com.hartwig.hmftools.lilac.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.registerCommonConfig;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.CLASS_1_EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.DEFLATE_TEMPLATE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.ReferenceData.loadHlaTranscripts;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_Y;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.EXON_BOUNDARY;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.IDENTICAL;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILDCARD;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.SEQUENCE_DELIM;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.buildAminoAcidSequenceFromNucleotides;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.lilac.CohortFrequency;
import com.hartwig.hmftools.lilac.GeneCache;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.lang3.NotImplementedException;
import org.jetbrains.annotations.NotNull;

public class GenerateReferenceSequences
{
    private static final String ALLELE_FREQUENCIES_FILENAME = "lilac_allele_frequencies.csv";
    private static final double WILDCARD_FREQUENCY_CUTOFF = 0.001;

    private final String mResourceDir;
    private final GeneCache mGeneCache;

    private final LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> mGroupedNucleotideSequences_;
    private final LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> mGroupedAminoAcidSequences_;

    private final HlaAlleleCache mAlleleCache;
    private final CohortFrequency mAlleleFrequencies_;

    public GenerateReferenceSequences(final ConfigBuilder configBuilder)
    {
        File alleleFrequenciesFile = new File(parseOutputDir(configBuilder), ALLELE_FREQUENCIES_FILENAME.toString());

        mAlleleCache = new HlaAlleleCache();
        mAlleleFrequencies_ = new CohortFrequency(alleleFrequenciesFile.toString());
        mResourceDir = configBuilder.getValue(RESOURCE_DIR);

        Map<HlaGene, TranscriptData> hlaTranscriptMap = loadHlaTranscripts(V37, null);
        mGeneCache = new GeneCache(hlaTranscriptMap);

        mGroupedNucleotideSequences_ = Maps.newLinkedHashMap();
        mGroupedAminoAcidSequences_ = Maps.newLinkedHashMap();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        registerCommonConfig(configBuilder);

        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateReferenceSequences seqGenerator = new GenerateReferenceSequences(configBuilder);

        LL_LOGGER.info("generating HLA sequence reference data");

        CLASS_1_EXCLUDED_ALLELES.clear();

        seqGenerator.loadSequenceFiles();

        final File outputDir = new File(parseOutputDir(configBuilder));
        seqGenerator.rewriteRefData(outputDir);

        LL_LOGGER.info("reference data written");
    }

    private void rewriteRefData(final File outputDir)
    {
        File nucleotideSequenceFile = new File(outputDir, NUC_REF_FILE);
        List<HlaSequenceLoci> nucleotideSequences = mGroupedNucleotideSequences_.values().stream()
                .flatMap(x -> x.stream())
                .toList();
        writeSequenceData(nucleotideSequenceFile, nucleotideSequences);

        File aminoAcidSequenceFile = new File(outputDir, AA_REF_FILE);
        List<HlaSequenceLoci> aminoAcidSequences = mGroupedAminoAcidSequences_.values().stream()
                .flatMap(x -> x.stream())
                .toList();
        writeSequenceData(aminoAcidSequenceFile, aminoAcidSequences);

        // TODO: GENE_CACHE is null
//        File aminoAcidComparisonFile = new File(outputDir, "hla_ref_aminoacid_compare.csv");
//        writeAminoAcidSequences(aminoAcidComparisonFile, mAminoAcidSequences);
    }

    private void addSequences(Map<HlaAllele, List<HlaSequenceLoci>> sequences, Map<HlaAllele, List<HlaSequenceLoci>> newSequences)
    {
        for(Map.Entry<HlaAllele, List<HlaSequenceLoci>> entry : newSequences.entrySet())
            sequences.computeIfAbsent(entry.getKey(), k -> Lists.newArrayList()).addAll(entry.getValue());
    }

    private void loadSequenceFiles()
    {
        LL_LOGGER.info("reading nucleotide files");

        for(HlaGene gene : mGeneCache.GeneNames)
            addSequences(mGroupedNucleotideSequences_, nucleotideLoci(new File(mResourceDir, gene.shortName() + "_nuc.txt"), gene));

        LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> yNucSequences = nucleotideLoci(new File(mResourceDir, "Y_nuc.txt"), HLA_Y);
        addSequences(mGroupedNucleotideSequences_, yNucSequences);

        // TODO: We do not include HLA-H seqs?
        // List<HlaSequenceLoci> hNucSequences = nucleotideLoci(mResourceDir + "H_nuc.txt");
        // mNucleotideSequences.addAll(hNucSequences);

        LL_LOGGER.info("reading protein files");

        for(HlaGene gene : mGeneCache.GeneNames)
            addSequences(mGroupedAminoAcidSequences_, aminoAcidLoci(new File(mResourceDir, gene.shortName() + "_prot.txt"), gene));

        HlaSequenceLoci sequenceTemplate = mGroupedAminoAcidSequences_.values().stream()
                .flatMap(x -> x.stream())
                .filter(x -> x.Allele.matches(DEFLATE_TEMPLATE))
                .findFirst()
                .orElse(null);

        // TODO: do we use H sequences?
        // build H and Y amino-acid sequences from their nucleotides
        if(sequenceTemplate != null)
        {
            yNucSequences.values().stream().flatMap(x -> x.stream()).forEach(x ->
            {
                HlaSequenceLoci aminoAcid = buildAminoAcidSequenceFromNucleotides(x, sequenceTemplate);
                mGroupedAminoAcidSequences_.computeIfAbsent(
                        aminoAcid.Allele.asAlleleGroup(), k -> Lists.newArrayList()).add(aminoAcid);
            });
            // hNucSequences.forEach(x -> mAminoAcidSequences.add(buildAminoAcidSequenceFromNucleotides(x, sequenceTemplate)));
        }
    }

    private LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> nucleotideLoci(final File filename, final HlaGene gene)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readDefintionFile(filename, gene);

        List<HlaSequence> filteredSequences = Lists.newArrayList(sequences);

        List<HlaSequence> reducedSequences = reduceToSixDigit(filteredSequences);

        List<HlaSequenceLoci> builtSequences = buildSequences(reducedSequences, false);
        LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> groupedBuiltSequences = Maps.newLinkedHashMap();
        for(HlaSequenceLoci seq : builtSequences)
            groupedBuiltSequences.computeIfAbsent(seq.Allele.asAlleleGroup(), k -> Lists.newArrayList()).add(seq);

        return fillWildcards(groupedBuiltSequences);
    }

    private LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> aminoAcidLoci(final File filename, final HlaGene gene)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readDefintionFile(filename, gene);

        List<HlaSequence> filteredSequences = Lists.newArrayList(sequences);

        List<HlaSequence> reducedSequences = reduceToFourDigit(filteredSequences);

        reducedSequences = reducedSequences.stream()
                .map(x -> x.getRawSequence().endsWith("X") ? x : x.copyWithAdditionalSequence("X"))
                .collect(Collectors.toList());

        List<HlaSequenceLoci> builtSequences = buildSequences(reducedSequences, true);
        LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> groupedBuiltSequences = Maps.newLinkedHashMap();
        for(HlaSequenceLoci seq : builtSequences)
            groupedBuiltSequences.computeIfAbsent(seq.Allele.asAlleleGroup(), k -> Lists.newArrayList()).add(seq);

        return fillWildcards(groupedBuiltSequences);
    }

    private LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> fillWildcards(final LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> sequences)
    {
        LinkedHashMap<HlaAllele, List<HlaSequenceLoci>> output = Maps.newLinkedHashMap();
        for(Map.Entry<HlaAllele, List<HlaSequenceLoci>> entry : sequences.entrySet())
        {
            HlaAllele alleleGroup = entry.getKey();
            List<HlaSequenceLoci> seqs = entry.getValue();

            // TODO:
            if(seqs.stream().anyMatch(x -> x.Allele.toString().equals("A*01:01:12")))
            {
                System.out.println("");
            }

            List<HlaSequenceLoci> filteredSeqs = Lists.newArrayList();
            double maxFrequency = -1.0;
            HlaSequenceLoci maxSeq = null;
            for(HlaSequenceLoci seq : seqs)
            {
                double frequency = mAlleleFrequencies_.getAlleleFrequency(seq.Allele.asFourDigit());
                if(frequency > maxFrequency)
                {
                    maxFrequency = frequency;
                    maxSeq = seq;
                }

                if(frequency < WILDCARD_FREQUENCY_CUTOFF)
                    continue;

                if(seq.hasWildcards())
                    continue;

                if(seq.Allele.asFourDigit().toString().matches("^.*[^0-9]$"))
                    continue;

                filteredSeqs.add(seq);
            }

            // TODO:
            if(maxSeq == null)
            {
                throw new IllegalStateException("maxSeq is null");
            }

            List<String> consensusSeq = filteredSeqs.isEmpty() ? maxSeq.getSequences() : buildConsensus(filteredSeqs);
            List<HlaSequenceLoci> newSeqs = seqs.stream().map(x -> fillWithConsensus(consensusSeq, x)).toList();
            output.put(alleleGroup, newSeqs);
        }

        return output;
    }

    private static List<String> buildConsensus(final List<HlaSequenceLoci> sequences)
    {
        // TODO:
        if(sequences.stream().anyMatch(x -> x.Allele.toString().equals("A*01:01:12")))
        {
            System.out.println("");
        }

        int maxLength = sequences.stream().mapToInt(x -> x.length()).max().getAsInt();

        // reverse the sequence building step
        final char sep = '@';
        List<StringBuilder> reversedSequences = Lists.newArrayList();
        for(int i = 0; i < sequences.size(); i++)
            reversedSequences.add(new StringBuilder());

        for(int i = 0; i < maxLength; i++)
        {
            final int index = i;
            List<String> loci = sequences.stream()
                    .map(x -> x.getSequences())
                    .map(x -> index >= x.size() ? null : x.get(index))
                    .toList();

            int maxLociLength = loci.stream().mapToInt(x -> x == null ? 0 : x.length()).max().getAsInt();
            for(int j = 0; j < loci.size(); j++)
            {
                String locus = loci.get(j);
                if(locus == null)
                    continue;

                reversedSequences.get(j).append(locus);
                if(locus.length() < maxLociLength)
                    reversedSequences.get(j).append(DEL_STR.repeat(maxLociLength - locus.length()));

                reversedSequences.get(j).append(sep);
            }
        }

        // now build consensus
        maxLength = reversedSequences.stream().mapToInt(x -> x.length()).max().getAsInt();
        StringBuilder consensusSeqBuilder = new StringBuilder();
        for(int i = 0; i < maxLength; i++)
        {
            Character consensusChar = null;
            boolean containsSep = false;
            for(int j = 0; j < reversedSequences.size(); j++)
            {
                StringBuilder seq = reversedSequences.get(j);
                if(i >= seq.length())
                    continue;

                // TODO:
                if(seq.charAt(i) == sep)
                    containsSep = true;

                if(consensusChar == null)
                {
                    consensusChar = seq.charAt(i);
                    continue;
                }

                if(!consensusChar.equals(seq.charAt(i)))
                {
                    consensusChar = WILDCARD;
                    break;
                }
            }

            // TODO:
            if(consensusChar == null)
            {
                throw new IllegalStateException("consensusChar is null");
            }

            // TODO:
            if(containsSep && !consensusChar.equals(sep))
            {
                throw new IllegalStateException("consensusChar is " + consensusChar + " but containsSep " + sep);
            }

            consensusSeqBuilder.append(consensusChar);
        }

        List<String> consensusSeqs = Lists.newArrayList();
        for(String seq : consensusSeqBuilder.toString().split("[" + sep + "]"))
        {
            if(seq.equals(DEL_STR))
            {
                consensusSeqs.add(seq);
                continue;
            }

            consensusSeqs.add(seq.replaceAll("[.]+$", ""));
        }

        return consensusSeqs;
    }

    private static HlaSequenceLoci fillWithConsensus(final List<String> consensus, final HlaSequenceLoci sequence)
    {
        // TODO:
        if(sequence.Allele.toString().equals("A*01:01:12"))
        {
            System.out.println("");
        }

        List<String> newSeqs = Lists.newArrayList();
        for(int i = 0; i < sequence.length(); i++)
        {
            if(i >= consensus.size())
            {
                newSeqs.add(sequence.sequence(i));
                continue;
            }

            if(sequence.sequence(i).equals(WILD_STR) && !consensus.get(i).equals(WILD_STR))
            {
                newSeqs.add(consensus.get(i));
                continue;
            }

            newSeqs.add(sequence.sequence(i));
        }

        return new HlaSequenceLoci(sequence.Allele, newSeqs);
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

    private static void writeSequenceData(final File filename_, final Iterable<HlaSequenceLoci> sequenceData)
    {
        try
        {
            final BufferedWriter writer = createBufferedWriter(filename_.toString(), false);

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
            LL_LOGGER.error("failed to write file {}: {}", filename_.toString(), e.toString());
        }
    }

    private void writeAminoAcidSequences(final File filename_, final Collection<HlaSequenceLoci> sequenceData)
    {
        // write each allele's amino acid sequence with reference to the template allele A*01:01

        HlaSequenceLoci refSequence = sequenceData.stream().filter(x -> x.Allele.matches(DEFLATE_TEMPLATE)).findFirst().orElse(null);

        if(refSequence == null)
        {
            LL_LOGGER.error("default template({}) not found");
            return;
        }

        List<Integer> sequenceMaxLengths = Lists.newArrayList();
        int sequenceLociMax = sequenceData.stream().mapToInt(HlaSequenceLoci::length).max().orElse(0);

        int maxAlleleNameLength = sequenceData.stream().mapToInt(x -> x.Allele.toString().length()).max().orElse(0);

        for(int i = 0; i < sequenceLociMax; ++i)
        {
            final int locus = i;

            int maxLocusLength = sequenceData.stream()
                    .filter(x -> locus < x.length())
                    .mapToInt(x -> x.sequence(locus).length()).max().orElse(0);

            sequenceMaxLengths.add(maxLocusLength);
        }

        int alleleNameLength = maxAlleleNameLength + 5;
        int maxInsertLen = 2;

        try
        {
            final BufferedWriter writer = createBufferedWriter(filename_.toString(), false);

            // write locus values first
            String titleStr = padString("Index", alleleNameLength);
            writer.write(titleStr);

            int segmentInsertCount = 0;

            writer.write(String.format("%03d", 0));

            for(int i = 1; i < sequenceMaxLengths.size(); ++i)
            {
                if(sequenceMaxLengths.get(i) > 1)
                    ++segmentInsertCount;

                if((i % 10) == 0)
                {
                    String padding = padString("", 7 + (maxInsertLen - 1) * segmentInsertCount);
                    writer.write(String.format("%s%03d", padding, i));
                    segmentInsertCount = 0;
                }
            }

            writer.newLine();
            writer.newLine();

            // then actual transcript positions for the exon boundaries
            for(HlaGene gene : mGeneCache.GeneNames)
            {
                TranscriptData transData = mGeneCache.GeneTranscriptMap.get(gene);
                String geneStr = padString(gene.shortName(), alleleNameLength);
                writer.write(geneStr);

                if(transData.Strand == POS_STRAND)
                {
                    for(int i = 0; i < transData.exons().size(); ++i)
                    {
                        ExonData exon = transData.exons().get(i);

                        if(transData.CodingStart > exon.End)
                            continue;
                        else if(transData.CodingEnd < exon.Start)
                            break;

                        writer.write(String.format("%d: %d - %d\t",
                                exon.Rank, max(exon.Start, transData.CodingStart), min(exon.End, transData.CodingEnd)));
                    }
                }
                else
                {
                    for(int i = transData.exons().size() - 1; i >= 0; --i)
                    {
                        ExonData exon = transData.exons().get(i);

                        if(transData.CodingEnd < exon.Start)
                            continue;
                        else if(transData.CodingStart > exon.End)
                            break;

                        writer.write(String.format("%d: %d - %d\t",
                                exon.Rank, max(exon.Start, transData.CodingStart), min(exon.End, transData.CodingEnd)));
                    }
                }

                writer.newLine();
            }

            writer.newLine();

            for(HlaSequenceLoci sequenceLoci : sequenceData)
            {
                HlaAllele allele = sequenceLoci.Allele;
                boolean isReference = sequenceLoci == refSequence;
                String alleleStr = padString(allele.toString(), alleleNameLength);
                writer.write(alleleStr);

                List<Integer> exonBoundaries = getAminoAcidExonBoundaries(allele.Gene);

                int nextExonIndex = 0;
                int nextExonBoundary = exonBoundaries.get(nextExonIndex);
                ++nextExonIndex;

                for(int locus = 0; locus < sequenceLoci.length(); ++locus)
                {
                    if(locus > nextExonBoundary)
                    {
                        writer.write(EXON_BOUNDARY);
                        if(nextExonIndex < exonBoundaries.size())
                        {
                            nextExonBoundary = exonBoundaries.get(nextExonIndex);
                            ++nextExonIndex;
                        }
                        else
                        {
                            nextExonBoundary = sequenceLociMax + 1;
                        }
                    }

                    String seq = sequenceLoci.sequence(locus);

                    if(!isReference)
                    {
                        String refSeq = locus < refSequence.length() ? refSequence.sequence(locus) : "";
                        if(refSeq.equals(seq))
                            seq = String.valueOf(IDENTICAL);
                    }

                    int maxSeqLength = min(sequenceMaxLengths.get(locus), maxInsertLen); // only show inserts with 1-extra base
                    seq = padString(seq, maxSeqLength);
                    writer.write(seq);
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write file {}: {}", filename_.toString(), e.toString());
        }
    }

    private static String padString(final String str, int reqLength)
    {
        return padString(str, reqLength, ' ');
    }

    private static String padString(final String str, int reqLength, char padChar)
    {
        if(str.length() >= reqLength)
            return str;

        return str + String.valueOf(padChar).repeat(reqLength - str.length());
    }
}
