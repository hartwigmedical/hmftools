package com.hartwig.hmftools.viridian.reference;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
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
    private final Map<String, ViralContig> mContigsByName;
    private final Map<OncologyGroup, OncologyGroupInfo> mOncologyGroupInfo;

    public VirusReference(
            List<ViralContig> contigs, SAMSequenceDictionary sequenceDictionary,
            Map<OncologyGroup, OncologyGroupInfo> oncologyGroupInfo)
    {
        mSequenceDictionary = sequenceDictionary;
        mContigsByName = contigs.stream().collect(toMap(ViralContig::name, identity()));
        mOncologyGroupInfo = oncologyGroupInfo;
    }

    public ViralContig contig(String name)
    {
        ViralContig contig = mContigsByName.get(name);
        if(contig == null)
        {
            throw new IllegalArgumentException("Unknown virus contig: " + name);
        }
        return contig;
    }

    // Contigs in FASTA order, matching the BWA index so an alignment's reference index resolves to a contig by position.
    public SAMSequenceDictionary sequenceDictionary()
    {
        return mSequenceDictionary;
    }

    public Set<OncologyGroup> oncologyGroups()
    {
        return mContigsByName.values().stream().map(ViralContig::oncologyGroup).collect(toSet());
    }

    public OncologyGroupInfo oncologyGroupInfo(OncologyGroup group)
    {
        OncologyGroupInfo info = mOncologyGroupInfo.get(group);
        if(info == null)
        {
            throw new IllegalArgumentException("Unknown oncology group: " + group);
        }
        return info;
    }

    public static VirusReference load(String fastaFile, String infoTsvFile, String oncologyGroupInfoTsvFile)
    {
        List<VirusInfo> info = VirusInfo.load(infoTsvFile);
        SAMSequenceDictionary dictionary = loadSequenceDictionary(fastaFile);
        List<ViralContig> contigs = joinFastaAndInfo(dictionary, info);
        Set<OncologyGroup> groups = contigs.stream().map(ViralContig::oncologyGroup).collect(toSet());
        Map<OncologyGroup, OncologyGroupInfo> oncologyGroupInfo = OncologyGroupInfo.load(oncologyGroupInfoTsvFile);
        if(!oncologyGroupInfo.keySet().equals(groups))
        {
            throw new UserInputError("Oncology group info groups do not match the virus reference groups");
        }
        return new VirusReference(contigs, dictionary, oncologyGroupInfo);
    }

    // Joins FASTA contigs to their info rows. Result in FASTA order.
    static List<ViralContig> joinFastaAndInfo(SAMSequenceDictionary dictionary, List<VirusInfo> info)
    {
        Map<String, VirusInfo> remainingInfo = new LinkedHashMap<>();
        for(VirusInfo virusInfo : info)
        {
            remainingInfo.put(virusInfo.contigName(), virusInfo);
        }
        List<ViralContig> result = new ArrayList<>();
        for(SAMSequenceRecord sequence : dictionary.getSequences())
        {
            String contig = sequence.getSequenceName();
            VirusInfo row = remainingInfo.remove(contig);
            if(row == null)
            {
                throw new UserInputError(String.format("Virus reference contig has no info row: %s", contig));
            }
            result.add(new ViralContig(contig, sequence.getSequenceLength(), row.virusName(), row.oncologyGroup()));
        }

        if(!remainingInfo.isEmpty())
        {
            throw new UserInputError(String.format("Virus reference info rows have no FASTA contig: %s", remainingInfo.keySet()));
        }

        return result;
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
