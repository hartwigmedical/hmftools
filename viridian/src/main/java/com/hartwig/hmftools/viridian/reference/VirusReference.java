package com.hartwig.hmftools.viridian.reference;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.common.UserInputError;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// The curated resource of virus genomes, virus information, and virus reporting policy data.
public class VirusReference
{
    private final SAMSequenceDictionary mSequenceDictionary;
    private final Map<String, VirusGenome> mGenomesByContigName;
    private final Map<VirusGroup, VirusGroupInfo> mVirusGroupInfo;

    public VirusReference(
            List<VirusGenome> genomes, SAMSequenceDictionary sequenceDictionary,
            Map<VirusGroup, VirusGroupInfo> virusGroupInfo)
    {
        mSequenceDictionary = sequenceDictionary;
        mGenomesByContigName = genomes.stream().collect(toMap(VirusGenome::contigName, identity()));
        mVirusGroupInfo = virusGroupInfo;
    }

    public VirusGenome genome(String contigName)
    {
        VirusGenome genome = mGenomesByContigName.get(contigName);
        if(genome == null)
        {
            throw new IllegalArgumentException("Unknown virus genome: " + contigName);
        }
        return genome;
    }

    // Contigs in FASTA order, matching the BWA index so an alignment's reference index resolves to a contig by position.
    public SAMSequenceDictionary sequenceDictionary()
    {
        return mSequenceDictionary;
    }

    public VirusGroupInfo virusGroupInfo(VirusGroup group)
    {
        VirusGroupInfo info = mVirusGroupInfo.get(group);
        if(info == null)
        {
            throw new IllegalArgumentException("Unknown virus group: " + group);
        }
        return info;
    }

    public static VirusReference load(String virusGenomesFasta, String virusGenomesTsv, String virusGroupsTsv)
    {
        List<VirusGenomeInfo> info = VirusGenomeInfo.load(virusGenomesTsv);
        SAMSequenceDictionary dictionary = loadSequenceDictionary(virusGenomesFasta);
        List<VirusGenome> genomes = joinFastaAndInfo(dictionary, info);
        Set<VirusGroup> groups = genomes.stream().map(VirusGenome::virusGroup).collect(toSet());
        Map<VirusGroup, VirusGroupInfo> virusGroupInfo = VirusGroupInfo.load(virusGroupsTsv);
        if(!virusGroupInfo.keySet().equals(groups))
        {
            throw new UserInputError(
                    "Virus groups in the virus groups TSV do not match those in the virus genomes TSV");
        }
        return new VirusReference(genomes, dictionary, virusGroupInfo);
    }

    // Joins FASTA contigs to their info rows. Result in FASTA order.
    static List<VirusGenome> joinFastaAndInfo(SAMSequenceDictionary dictionary, List<VirusGenomeInfo> info)
    {
        Map<String, Virus> virusesByContigName = info.stream()
                .collect(toMap(VirusGenomeInfo::contigName, VirusGenomeInfo::virus));
        Set<String> fastaContigs = dictionary.getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName)
                .collect(toSet());

        List<String> contigsWithoutInfo = fastaContigs.stream()
                .filter(contig -> !virusesByContigName.containsKey(contig)).sorted().toList();
        List<String> infoWithoutContig = virusesByContigName.keySet().stream()
                .filter(contig -> !fastaContigs.contains(contig)).sorted().toList();
        if(!contigsWithoutInfo.isEmpty() || !infoWithoutContig.isEmpty())
        {
            throw new UserInputError(String.format(
                    "Virus reference FASTA contigs and info rows differ. Contigs without info row: %s. Info rows without contig: %s",
                    contigsWithoutInfo, infoWithoutContig));
        }

        return dictionary.getSequences().stream()
                .map(sequence -> new VirusGenome(
                        sequence.getSequenceName(), sequence.getSequenceLength(),
                        virusesByContigName.get(sequence.getSequenceName())))
                .toList();
    }

    private static SAMSequenceDictionary loadSequenceDictionary(String fastaFile)
    {
        try(IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(fastaFile)))
        {
            SAMSequenceDictionary dictionary = fasta.getSequenceDictionary();
            if(dictionary == null)
            {
                throw new UserInputError("Virus reference FASTA has no sequence dictionary (.dict): " + fastaFile);
            }
            return dictionary;
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to read virus reference FASTA index", e);
        }
    }

}
