package com.hartwig.hmftools.common.cider;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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
        primaryAssembly,
        chromosome,
        posStart,
        posEnd,
        strand,
        anchorStart,
        anchorEnd,
        anchorSequence,
        anchorAA
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
                boolean isPrimaryAssembly = record.getBoolean(Column.primaryAssembly);

                String anchorSequence = record.get(Column.anchorSequence);
                if(anchorSequence.isEmpty())
                {
                    anchorSequence = null;
                }

                String chromosome = record.get(Column.chromosome).intern();
                ChrBaseRegion genomicLocation = null;
                ChrBaseRegion anchorLocation = null;
                Strand strand = null;

                if(!chromosome.isEmpty())
                {
                    chromosome = refGenomeVersion.versionedChromosome(chromosome);

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

                    genomicLocation = new ChrBaseRegion(chromosome, posStart, posEnd);

                    String anchorStart = record.get(Column.anchorStart);
                    String anchorEnd = record.get(Column.anchorEnd);

                    if(!anchorStart.isEmpty() && !anchorEnd.isEmpty())
                    {
                        anchorLocation = new ChrBaseRegion(chromosome, Integer.parseInt(anchorStart), Integer.parseInt(anchorEnd));
                    }
                }

                igTcrGenes.add(new IgTcrGene(geneName, allele, region, functionality, genomicLocation, strand, isPrimaryAssembly,
                        anchorSequence, anchorLocation));
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
                    row.set(Column.primaryAssembly, gene.inPrimaryAssembly());
                    row.set(Column.chromosome, gene.geneLocation() != null ? gene.geneLocation().chromosome() : "");
                    row.set(Column.posStart, gene.geneLocation() != null ? Integer.toString(gene.geneLocation().start()) : "");
                    row.set(Column.posEnd, gene.geneLocation() != null ? Integer.toString(gene.geneLocation().end()) : "");
                    row.set(Column.strand, gene.geneStrand() != null ? String.valueOf(gene.geneStrand().asChar()) : "");
                    row.set(Column.anchorStart, gene.anchorLocation() != null ? Integer.toString(gene.anchorLocation().start()) : "");
                    row.set(Column.anchorEnd, gene.anchorLocation() != null ? Integer.toString(gene.anchorLocation().end()) : "");
                    row.set(Column.anchorSequence, gene.anchorSequence() != null ? gene.anchorSequence() : "");
                    row.set(Column.anchorAA, gene.anchorSequence() != null ? Codons.aminoAcidFromBases(gene.anchorSequence()) : "");
                }
        );
    }
}
