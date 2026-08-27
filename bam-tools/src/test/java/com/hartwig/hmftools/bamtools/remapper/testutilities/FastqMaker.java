package com.hartwig.hmftools.bamtools.remapper.testutilities;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class FastqMaker
{
    public static void main(String[] args) throws Exception
    {
        FastqMaker maker = new FastqMaker(HumanChromosome._7, 44000000, 45000000);
        final File destination = new File("/Users/timlavers/work/batches/2026/7/17/2/chr7_44_45.fastq");
        Preconditions.checkArgument(destination.getParentFile().exists());
        Preconditions.checkArgument(!destination.exists());
        maker.writeFasta(destination);
    }

    private final HumanChromosome chromosome;
    private final int start;
    private final int end;
    private final String qualityString = "F".repeat(150);

    public FastqMaker(final HumanChromosome chromosome, final int start, final int end)
    {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    public void writeFasta(File destination) throws Exception
    {
        try(var writer = Files.newBufferedWriter(destination.toPath(), StandardCharsets.UTF_8))
        {
            final String chrName = V38.versionedChromosome(chromosome);
            final File refGenomeFile = new File("/Users/timlavers/work/data/reference_genomes/38/Homo_sapiens_assembly38.alt.masked.fasta");
            RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
            int currentStart = start;
            int currentEnd = currentStart + 450;
            while(currentEnd < end)
            {
                var chunk = refGenomeSource.getBaseString(chrName, currentStart, currentEnd);
                var readNamePrefix = "@" + chrName + "_" + currentStart + "_" + currentEnd;
                writer.write(readNamePrefix + " 1\n");
                var firstRead = chunk.substring(0, 150);
                writer.write(firstRead + "\n");
                writer.write("+\n");
                writer.write(qualityString + "\n");

                var secondRead = Nucleotides.reverseComplementBases(chunk.substring(300, 450));
                writer.write(readNamePrefix + " 2\n");
                writer.write(secondRead + "\n");
                writer.write("+\n");
                writer.write(qualityString + "\n");

                currentStart += 150;
                currentEnd += 150;
            }
        }
    }
}
