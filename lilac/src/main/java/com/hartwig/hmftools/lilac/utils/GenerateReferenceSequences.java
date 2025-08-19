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
import static com.hartwig.hmftools.lilac.LilacConfig.MHC_CLASS;
import static com.hartwig.hmftools.lilac.LilacConfig.RESOURCE_DIR;
import static com.hartwig.hmftools.lilac.LilacConfig.registerCommonConfig;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.CLASS_1_EXCLUDED_ALLELES;
import static com.hartwig.hmftools.lilac.ReferenceData.AA_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.DEFLATE_TEMPLATE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_REF_FILE;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.ReferenceData.loadHlaTranscripts;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.EXON_BOUNDARY;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.IDENTICAL;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.SEQUENCE_DELIM;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.buildAminoAcidSequenceFromNucleotides;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.lilac.GeneCache;
import com.hartwig.hmftools.lilac.MhcClass;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaAlleleCache;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.jetbrains.annotations.NotNull;

public class GenerateReferenceSequences
{
    private final String mResourceDir;
    private final GeneCache mGeneCache;

    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mAminoAcidSequences;

    private final HlaAlleleCache mAlleleCache;

    public GenerateReferenceSequences(final ConfigBuilder configBuilder)
    {
        mAlleleCache = new HlaAlleleCache();
        mResourceDir = configBuilder.getValue(RESOURCE_DIR);

        MhcClass classType = MhcClass.valueOf(configBuilder.getValue(MHC_CLASS));

        if(classType != MhcClass.CLASS_1)
        {
            LL_LOGGER.error("only class-1 support for now");
            System.exit(1);
        }

        Map<HlaGene, TranscriptData> hlaTranscriptMap = loadHlaTranscripts(V37, classType);
        mGeneCache = new GeneCache(classType, hlaTranscriptMap);

        mNucleotideSequences = Lists.newArrayList();
        mAminoAcidSequences = Lists.newArrayList();
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

        if(!seqGenerator.loadSequenceFiles())
        {
            LL_LOGGER.error("ref data loading failed");
            System.exit(1);
        }

        final String outputDir = parseOutputDir(configBuilder);
        seqGenerator.rewriteRefData(outputDir);

        LL_LOGGER.info("reference data written");
    }

    private void rewriteRefData(final String outputDir)
    {
        String nucleotideSequenceFile = outputDir + NUC_REF_FILE;
        writeSequenceData(nucleotideSequenceFile, mNucleotideSequences);

        String aminoAcidSequenceFile = outputDir + AA_REF_FILE;
        writeSequenceData(aminoAcidSequenceFile, mAminoAcidSequences);

        String aminoAcidComparisonFile = outputDir + "hla_ref_aminoacid_compare.csv";
        writeAminoAcidSequences(aminoAcidComparisonFile, mAminoAcidSequences);
    }

    private boolean loadSequenceFiles()
    {
        LL_LOGGER.info("reading nucleotide files");

        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "A_nuc.txt"));
        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "B_nuc.txt"));
        mNucleotideSequences.addAll(nucleotideLoci(mResourceDir + "C_nuc.txt"));

        List<HlaSequenceLoci> yNucSequences = nucleotideLoci(mResourceDir + "Y_nuc.txt");
        mNucleotideSequences.addAll(yNucSequences);

        // List<HlaSequenceLoci> hNucSequences = nucleotideLoci(mResourceDir + "H_nuc.txt");
        // mNucleotideSequences.addAll(hNucSequences);

        LL_LOGGER.info("reading protein files");

        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "A_prot.txt"));
        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "B_prot.txt"));
        mAminoAcidSequences.addAll(aminoAcidLoci(mResourceDir + "C_prot.txt"));

        HlaSequenceLoci sequenceTemplate = mAminoAcidSequences.stream()
                .filter(x -> x.Allele.matches(DEFLATE_TEMPLATE)).findFirst().orElse(null);

        // build H and Y amino-acid sequences from their nucleotides
        if(sequenceTemplate != null)
        {
            yNucSequences.forEach(x -> mAminoAcidSequences.add(buildAminoAcidSequenceFromNucleotides(x, sequenceTemplate)));
            // hNucSequences.forEach(x -> mAminoAcidSequences.add(buildAminoAcidSequenceFromNucleotides(x, sequenceTemplate)));
        }

        return true;
    }

    private List<HlaSequenceLoci> nucleotideLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readDefintionFile(filename);

        List<HlaSequence> filteredSequences = Lists.newArrayList(sequences);

        List<HlaSequence> reducedSequences = reduceToSixDigit(filteredSequences);

        return buildSequences(reducedSequences, false);
    }

    private List<HlaSequenceLoci> aminoAcidLoci(final String filename)
    {
        List<HlaSequence> sequences = HlaSequenceFile.readDefintionFile(filename);

        List<HlaSequence> filteredSequences = Lists.newArrayList(sequences);

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

    private static void writeSequenceData(final String filename, final Iterable<HlaSequenceLoci> sequenceData)
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

    private void writeAminoAcidSequences(final String filename, final Collection<HlaSequenceLoci> sequenceData)
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
            final BufferedWriter writer = createBufferedWriter(filename, false);

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
            LL_LOGGER.error("failed to write file {}: {}", filename, e.toString());
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
