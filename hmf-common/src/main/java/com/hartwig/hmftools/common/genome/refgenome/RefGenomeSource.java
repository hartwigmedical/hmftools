package com.hartwig.hmftools.common.genome.refgenome;

import java.io.File;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeSource implements RefGenomeInterface
{
    private final IndexedFastaSequenceFile mRefGenome;

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeSource.class);

    public RefGenomeSource(final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
    }

    @Override
    public String getBaseString(final String chromosome, int posStart, int posEnd)
    {
        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBaseString();
    }

    public static RefGenomeSource loadRefGenome(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return null;

        try
        {
            LOGGER.debug("loading indexed fasta reference file");
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
