package com.hartwig.hmftools.cdr3;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.io.Files;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class Cdr3GeneLoader implements VJGeneStore
{
    private static final Logger sLogger = LogManager.getLogger(Cdr3GeneLoader.class);

    private final List<VJGene> mVJGenes;
    private final Multimap<String, VJGene> mAnchorSequenceMap = ArrayListMultimap.create();

    private final Map<VJGeneType, Multimap<String, VJGene>> mGeneTypeAnchorSeqMap = new HashMap<>();

    private final Multimap<VJAnchorReferenceLocation, VJGene> mGeneLocationVJGeneMap = ArrayListMultimap.create();

    @Override
    public List<VJGene> getVJGenes()
    {
        return mVJGenes;
    }

    @Override
    public Set<String> getAnchorSequenceSet()
    {
        return mAnchorSequenceMap.keySet();
    }

    @Override
    public Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType)
    {
        Multimap<String, VJGene> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.keySet() : Collections.emptySet();
    }

    @Override
    public Collection<VJGene> getByAnchorSequence(@NotNull String anchorSeq)
    {
        return mAnchorSequenceMap.get(anchorSeq);
    }

    @Override
    public Collection<VJGene> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq)
    {
        Multimap<String, VJGene> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.get(anchorSeq) : Collections.emptySet();
    }

    @Override
    public Collection<VJGene> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation)
    {
        return mGeneLocationVJGeneMap.get(vjAnchorReferenceLocation);
    }

    @Override
    public Collection<VJAnchorReferenceLocation> getVJAnchorReferenceLocations()
    {
        return mGeneLocationVJGeneMap.keySet();
    }

    public Cdr3GeneLoader(RefGenomeVersion refGenomeVersion) throws IOException
    {
        mVJGenes = loadImgtGeneTsv(refGenomeVersion);

        // from this we find all the anchor sequence locations and fix them
        for (VJGene gene : mVJGenes)
        {
            if (gene.getAnchorLocation() != null)
            {
                mGeneLocationVJGeneMap.put(new VJAnchorReferenceLocation(gene.getType().getVj(), gene.getAnchorLocation()), gene);
            }

            if (!gene.getAnchorSequence().isEmpty())
            {
                mAnchorSequenceMap.put(gene.getAnchorSequence(), gene);
                mGeneTypeAnchorSeqMap.computeIfAbsent(gene.getType(), o -> ArrayListMultimap.create())
                    .put(gene.getAnchorSequence(), gene);
            }
        }

        sLogger.info("found {} gene locations", mGeneLocationVJGeneMap.keySet().size());
    }

    public static List<ReferenceSequence> loadFasta(String resourcePath) throws IOException
    {
        List<ReferenceSequence> seqList = new ArrayList<>();
        try (java.io.InputStream fastaStream = Cdr3GeneLoader.class.getClassLoader().getResourceAsStream(resourcePath))
        {
            if (fastaStream == null)
            {
                sLogger.error("unable to find resource file: {}", resourcePath);
                throw new RuntimeException("unable to find resource file: " + resourcePath);
            }

            File tempFile = null;
            // write the resource out to a temp file and read it back
            try
            {
                tempFile = java.io.File.createTempFile(Files.getNameWithoutExtension(resourcePath), ".fa");
                Files.write(fastaStream.readAllBytes(), tempFile);
                try (var fastaFile = new FastaSequenceFile(tempFile, false))
                {
                    ReferenceSequence seq = fastaFile.nextSequence();
                    while (seq != null)
                    {
                        seqList.add(seq);
                        seq = fastaFile.nextSequence();
                    }
                }
            }
            finally
            {
                if (tempFile != null)
                {
                    //noinspection ResultOfMethodCallIgnored
                    tempFile.delete();
                }
            }
        }

        return seqList;
    }

    private static List<VJGene> loadImgtGeneTsv(RefGenomeVersion refGenomeVersion) throws IOException
    {
        String resourcePath;

        if (refGenomeVersion.is37())
        {
            resourcePath = "imgt_anchor.37.tsv";
        }
        else if (refGenomeVersion.is38())
        {
            resourcePath = "imgt_anchor.38.tsv";
        }
        else
        {
            throw new IllegalArgumentException("unknown ref genome version: " + refGenomeVersion.toString());
        }
        List<VJGene> VJGeneList = new ArrayList<>();

        java.io.InputStream tsvStream = Cdr3GeneLoader.class.getClassLoader().getResourceAsStream(resourcePath);
        if (tsvStream == null)
        {
            sLogger.error("unable to find resource file: {}", resourcePath);
            throw new RuntimeException("unable to find resource file: " + resourcePath);
        }

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(tsvStream)))
        {
            CSVFormat format = CSVFormat.Builder.create()
                .setDelimiter('\t')
                .setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                .build();
            Iterable<CSVRecord> records = format.parse(reader);
    
            for (CSVRecord record : records)
            {
                String id = record.get("id");
                String name = record.get("gene");
                String allele = record.get("allele");
                int posStart = Integer.parseInt(record.get("posStart"));
                int posEnd = Integer.parseInt(record.get("posEnd"));
                String anchorSequence = record.get("anchorSequence");
                String chromosome = record.get("chr");

                if (refGenomeVersion.is37())
                    chromosome = chromosome.replace("chr", "");

                String strandStr = record.get("strand");
                Strand strand = null;
                if (strandStr.equals("+"))
                    strand = Strand.FORWARD;
                else if (strandStr.equals("-"))
                    strand = Strand.REVERSE;
                int anchorStart = Integer.parseInt(record.get("anchorStart"));
                int anchorEnd = Integer.parseInt(record.get("anchorEnd"));
                String sequence = record.get("sequence");

                GeneLocation geneLocation = null;
                GeneLocation anchorLocation = null;

                if (!chromosome.isEmpty())
                {
                    if (posStart <= 0 || posEnd <= 0)
                    {
                        throw new RuntimeException("chromosome exist but pos start or pos end invalid");
                    }

                    if (strand == null)
                    {
                        throw new RuntimeException("chromosome exist but strand invalid");
                    }

                    geneLocation = new GeneLocation(chromosome, posStart, posEnd, strand);

                    if (anchorStart >= 0 && anchorEnd >= 0)
                        anchorLocation = new GeneLocation(chromosome, anchorStart, anchorEnd, strand);
                }

                VJGene VJGene = new VJGene(
                        id, name, allele, geneLocation, sequence, anchorSequence, anchorLocation);
                VJGeneList.add(VJGene);
            }
        }
        return VJGeneList;
    }
}
