package com.hartwig.hmftools.amber.e2e;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.bam.testutilities.PairedRecordsBuilder;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;

public class AmberReadsSpecification
{
    private final ListMultimap<Chromosome, AmberSiteRead> Sites = ArrayListMultimap.create();

    public AmberReadsSpecification(final ListMultimap<Chromosome, AmberSite> amberSites, final File readsFile) throws Exception
    {
        List<String> lines = FileUtils.readLines(readsFile, StandardCharsets.UTF_8);
        lines.forEach(line ->
        {
            // For now assume that there's just one site per line.
            String[] split = line.split(" ");
            AmberSiteReadsSpecification spec = new AmberSiteReadsSpecification(split[0].trim(), split[1].trim());
            Sites.putAll(spec.chromosome(), spec.resolve(amberSites));
        });
    }

    public void writeBam(RefGenomeSource refGenomeSource, String outputFileName) throws Exception
    {
        var header = new SAMFileHeader();
        Sites.keySet().forEach(chromosomeIndex ->
        {
            String chrName = RefGenomeVersion.V38.versionedChromosome(chromosomeIndex);
            int chrLength = refGenomeSource.getChromosomeLength(chrName);
            header.getSequenceDictionary().addSequence(new SAMSequenceRecord(chrName, chrLength));
        });
        SAMFileWriter bamWriter = new FastBamWriter(header, outputFileName);
        writeBamEntries(bamWriter, refGenomeSource);
        bamWriter.close();
    }

    private void writeBamEntries(SAMFileWriter bamWriter, RefGenomeSource refGenomeSource)
    {
        int readNumber = 0;
        for(Chromosome chromosome : Sites.keySet())
        {
            for(AmberSiteRead amberSiteRead : Sites.get(chromosome))
            {
                var readName = "A:B:C:" + readNumber;
                var bases = amberSiteRead.baseRegion(refGenomeSource);
                var mateBases = amberSiteRead.mateReadBasesRegion(refGenomeSource);
                var records = new PairedRecordsBuilder(readName, bamWriter.getFileHeader()).build(Pair.of(bases, mateBases));
                bamWriter.addAlignment(records.getRight());
                bamWriter.addAlignment(records.getLeft());
            }
        }
    }
}
