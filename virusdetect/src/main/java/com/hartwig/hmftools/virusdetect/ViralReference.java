package com.hartwig.hmftools.virusdetect;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// The curated viral reference: contig names and lengths (from the FASTA sequence dictionary) joined
// to each contig's virus name and oncology group (from the info TSV). Contigs are held in FASTA
// order, the alignment-target order relied on when translating BWA reference indices to names.
public class ViralReference
{
    private final List<ViralContig> mContigs;
    private final SAMSequenceDictionary mSequenceDictionary;

    private static final Logger LOGGER = LogManager.getLogger(ViralReference.class);

    private ViralReference(List<ViralContig> contigs, SAMSequenceDictionary sequenceDictionary)
    {
        mContigs = contigs;
        mSequenceDictionary = sequenceDictionary;
    }

    public List<ViralContig> contigs()
    {
        return mContigs;
    }

    // The FASTA sequence dictionary, in FASTA order, matching the BWA index: an alignment's reference
    // index resolves to a contig by position, and this is the header for the aligned BAM.
    public SAMSequenceDictionary sequenceDictionary()
    {
        return mSequenceDictionary;
    }

    public static ViralReference load(String fastaFile, String infoTsvFile)
    {
        Map<String, InfoRow> info = loadInfo(infoTsvFile);
        SAMSequenceDictionary dictionary = loadSequenceDictionary(fastaFile);
        List<ViralContig> contigs = join(dictionary, info);
        LOGGER.info("loaded viral reference: {} contigs", contigs.size());
        return new ViralReference(contigs, dictionary);
    }

    // Joins FASTA contigs (in FASTA order) to their info rows, enforcing the 1:1 integrity rule.
    static List<ViralContig> join(SAMSequenceDictionary dictionary, Map<String, InfoRow> info)
    {
        Map<String, InfoRow> remaining = new LinkedHashMap<>(info);
        List<ViralContig> contigs = new ArrayList<>();
        for(SAMSequenceRecord sequence : dictionary.getSequences())
        {
            String contig = sequence.getSequenceName();
            InfoRow row = remaining.remove(contig);
            if(row == null)
            {
                throw new UserInputError(String.format("viral reference contig has no info row: %s", contig));
            }
            contigs.add(new ViralContig(contig, sequence.getSequenceLength(), row.virusName(), row.oncologyGroup()));
        }

        if(!remaining.isEmpty())
        {
            throw new UserInputError(String.format("viral reference info rows have no FASTA contig: %s", remaining.keySet()));
        }

        return contigs;
    }

    static Map<String, InfoRow> loadInfo(String infoTsvFile)
    {
        Map<String, InfoRow> info = new LinkedHashMap<>();
        try(DelimFileReader reader = new DelimFileReader(infoTsvFile))
        {
            List<String> columns = reader.getColumnNames();
            for(Column column : Column.values())
            {
                if(!columns.contains(column.name()))
                {
                    throw new UserInputError(String.format("viral reference info missing column: %s", column.name()));
                }
            }

            for(DelimFileReader.Row row : reader)
            {
                String contig = row.get(Column.ref_contig);
                InfoRow previous = info.put(contig, new InfoRow(row.get(Column.virus_name), row.get(Column.oncology_group)));
                if(previous != null)
                {
                    throw new UserInputError(String.format("viral reference info has duplicate contig: %s", contig));
                }
            }
        }
        return info;
    }

    private static SAMSequenceDictionary loadSequenceDictionary(String fastaFile)
    {
        try(IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(fastaFile)))
        {
            SAMSequenceDictionary dictionary = fasta.getSequenceDictionary();
            if(dictionary == null)
            {
                throw new UserInputError("viral reference FASTA has no sequence dictionary (.dict): " + fastaFile);
            }
            return dictionary;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read viral reference FASTA index", e);
        }
    }

    private enum Column
    {
        ref_contig,
        virus_name,
        oncology_group
    }

    record InfoRow(String virusName, String oncologyGroup)
    {
    }
}
