package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;
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
        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBases();
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
}
