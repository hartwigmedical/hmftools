package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeSource implements RefGenomeInterface
{
    public static final String REF_GENOME = "ref_genome";
    public static final String REF_GENOME_CFG_DESC = "Path to reference genome fasta files";

    private final IndexedFastaSequenceFile mRefGenome;

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeSource.class);

    public IndexedFastaSequenceFile refGenomeFile() { return mRefGenome; }

    public static void addRefGenomeVersion(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
    }

    public static void addRefGenomeFile(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPath(REF_GENOME, required, REF_GENOME_CFG_DESC);
    }

    public static void addRefGenomeConfig(final ConfigBuilder configBuilder, boolean required)
    {
        if(required)
            configBuilder.addRequiredConfigItem(REF_GENOME_VERSION, REF_GENOME_VERSION_CFG_DESC);
        else
            configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());

        configBuilder.addPath(REF_GENOME, required, REF_GENOME_CFG_DESC);
    }

    public RefGenomeSource(final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
    }

    @Override
    public String getBaseString(final String chromosome, int posStart, int posEnd)
    {
        if(posStart < 1 || posEnd < posStart)
            return null;

        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBaseString();
    }

    @Override
    public String getBaseString(final String chromosome, final List<int[]> baseRanges)
    {
        StringBuilder refBases = new StringBuilder();
        baseRanges.forEach(x -> refBases.append(mRefGenome.getSubsequenceAt(chromosome, x[0], x[1]).getBaseString()));
        return refBases.toString();
    }

    @Override
    public int getChromosomeLength(final String chromosome)
    {
        return mRefGenome.getSequenceDictionary().getSequence(chromosome).getSequenceLength();
    }

    @Override
    public byte[] getBases(final String chromosome, int posStart, int posEnd)
    {
        if(posStart < 1 || posEnd < posStart)
            return null;

        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBases();
    }

    @Override
    public Map<String,Integer> chromosomeLengths()
    {
        Map<String,Integer> chromosomeLengthMap = Maps.newHashMap();

        for(SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            chromosomeLengthMap.put(sequenceRecord.getSequenceName(), sequenceRecord.getSequenceLength());
        }

        return chromosomeLengthMap;
    }

    public static RefGenomeVersion deriveRefGenomeVersion(final RefGenomeSource refGenomeSource)
    {
        String firstChromosome = refGenomeSource.refGenomeFile().getSequenceDictionary().getSequences().get(0).getSequenceName();
        return firstChromosome.startsWith(CHR_PREFIX) ? V38 : V37;
    }

    public static RefGenomeSource loadRefGenome(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return null;

        try
        {
            // LOGGER.debug("loading indexed fasta reference file");
            IndexedFastaSequenceFile refFastaSeqFile = new IndexedFastaSequenceFile(new File(filename));
            return new RefGenomeSource(refFastaSeqFile);
        }
        catch (IOException e)
        {
            LOGGER.error("Reference file loading failed: {}", e.toString());
            return null;
        }
    }

    public List<ChrBaseRegion> formRefGenomeRegions(final SpecificRegions specificRegions)
    {
        List<ChrBaseRegion> inputRegions = Lists.newArrayList();

        // form regions in same order as in dictionary
        for(SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            String contig = sequenceRecord.getSequenceName();

            if(specificRegions.excludeChromosome(contig))
                continue;

            inputRegions.add(new ChrBaseRegion(contig, 1, sequenceRecord.getSequenceLength()));
        }

        return inputRegions;
    }
}
