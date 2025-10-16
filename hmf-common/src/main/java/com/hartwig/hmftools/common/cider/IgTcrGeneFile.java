package com.hartwig.hmftools.common.cider;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class IgTcrGeneFile
{
    public enum Column
    {
        gene,
        allele,
        region,
        functionality,
        sequence,
        contig,
        posStart,
        posEnd,
        strand,
        anchorStart,
        anchorEnd,
        anchorSequence
    }

    private static final Logger sLogger = LogManager.getLogger(IgTcrGeneFile.class);

    public static List<IgTcrGene> read(RefGenomeVersion refGenomeVersion)
    {
        String resourcePath = "/refgenome/";
        if(refGenomeVersion.is37())
        {
            resourcePath += "igtcr_gene.37.tsv";
        }
        else if(refGenomeVersion.is38())
        {
            resourcePath += "igtcr_gene.38.tsv";
        }
        else
        {
            sLogger.error("unknown ref genome version: {}", refGenomeVersion);
            throw new IllegalArgumentException("unknown ref genome version: " + refGenomeVersion);
        }

        List<IgTcrGene> igTcrGenes = new ArrayList<>();

        InputStream tsvStream = IgTcrGeneFile.class.getResourceAsStream(resourcePath);
        if(tsvStream == null)
        {
            sLogger.error("unable to find resource file: {}", resourcePath);
            throw new RuntimeException("unable to find resource file: " + resourcePath);
        }

        try(DelimFileReader reader = new DelimFileReader(new BufferedReader(new InputStreamReader(tsvStream))))
        {
            for(DelimFileReader.Row record : reader)
            {
                String geneName = record.get(Column.gene);
                String allele = record.get(Column.allele);
                IgTcrRegion region = IgTcrRegion.valueOf(record.get(Column.region));
                IgTcrFunctionality functionality = IgTcrFunctionality.fromCode(record.get(Column.functionality));
                String sequence = record.get(Column.sequence);
                String anchorSequence = record.getOrNull(Column.anchorSequence);

                String contig = record.getOrNull(Column.contig);
                BaseRegion geneLocation = null;
                BaseRegion anchorLocation = null;
                Strand strand = null;

                if(contig != null)
                {
                    contig = refGenomeVersion.versionedChromosome(contig);

                    int posStart = record.getInt(Column.posStart);
                    int posEnd = record.getInt(Column.posEnd);

                    String strandStr = record.get(Column.strand);

                    if(strandStr.isEmpty())
                    {
                        throw new RuntimeException("chromosome exists but strand is missing");
                    }

                    strand = Strand.fromChar(strandStr.charAt(0));

                    if(posStart <= 0 || posEnd <= 0)
                    {
                        throw new RuntimeException("chromosome exists but pos start or pos end is invalid");
                    }

                    geneLocation = new BaseRegion(posStart, posEnd);

                    Integer anchorStart = record.getIntOrNull(Column.anchorStart);
                    Integer anchorEnd = record.getIntOrNull(Column.anchorEnd);

                    if(anchorStart != null && anchorEnd != null)
                    {
                        anchorLocation = new BaseRegion(anchorStart, anchorEnd);
                    }
                }

                igTcrGenes.add(new IgTcrGene(
                        geneName, allele, region, functionality, sequence, contig, geneLocation, strand, anchorSequence, anchorLocation));
            }
        }

        return igTcrGenes;
    }

    public static void write(String tsvPath, List<IgTcrGene> geneList)
    {
        DelimFileWriter.write(tsvPath, Column.values(), geneList,
                (gene, row) ->
                {
                    row.set(Column.gene, gene.geneName());
                    row.set(Column.allele, gene.allele());
                    row.set(Column.region, String.valueOf(gene.region()));
                    row.set(Column.functionality, gene.functionality().toCode());
                    row.set(Column.sequence, gene.sequence());
                    row.setOrNull(Column.contig, gene.contigName());
                    row.setOrNull(Column.posStart, gene.genePosition() != null ? gene.genePosition().start() : null);
                    row.setOrNull(Column.posEnd, gene.genePosition() != null ? gene.genePosition().end() : null);
                    row.setOrNull(Column.strand, gene.geneStrand() != null ? gene.geneStrand().asChar() : null);
                    row.setOrNull(Column.anchorStart, gene.anchorPosition() != null ? gene.anchorPosition().start() : null);
                    row.setOrNull(Column.anchorEnd, gene.anchorPosition() != null ? gene.anchorPosition().end() : null);
                    row.setOrNull(Column.anchorSequence, gene.anchorSequence());
                }
        );
    }
}
