package com.hartwig.hmftools.common.bamops;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class BamSampler
{
    private final BamSlicer mSlicer;
    private final int mMaxReadCount;
    private final RefGenomeSource mRefGenome;

    private int mReadCount;
    private int mMaxReadLength;
    private boolean mReadsPaired;
    private boolean mMateCigarSet;
    private Consumer<SAMRecord> mConsumer;

    private static final int DEFAULT_MAX_READS = 1000;

    public BamSampler(final String referenceGenome)
    {
        this(referenceGenome, DEFAULT_MAX_READS);
    }

    public BamSampler(final String referenceGenome, final int maxReadCount)
    {
        mReadCount = 0;
        mMaxReadLength = 0;
        mMaxReadCount = maxReadCount;
        mRefGenome = loadRefGenome(referenceGenome);
        mConsumer = null;

        mMateCigarSet = false;
        mReadsPaired = false;

        mSlicer = new BamSlicer(0);
    }

    public void setConsumer(final Consumer<SAMRecord> consumer) { mConsumer = consumer; }

    public int maxReadLength() { return mMaxReadLength; }
    public boolean readsPaired() { return mReadsPaired; }
    public boolean hasMateCigarSet() { return mMateCigarSet; }

    public ChrBaseRegion defaultRegion()
    {
        String sampleChromosome = mRefGenome.refGenomeFile().getSequenceDictionary().getSequence(0).getSequenceName();
        return new ChrBaseRegion(sampleChromosome, 100000, 200000);
    }

    public boolean calcBamCharacteristics(final String bamFile, final ChrBaseRegion sampleRegion)
    {
        if(mRefGenome == null)
            return false;

        if(!Files.exists(Paths.get(bamFile)))
            return false;

        SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(mRefGenome.refGenomeFile()))
                .open(new File(bamFile));

        // reset
        mReadCount = 0;
        mMaxReadLength = 0;

        mSlicer.slice(samReader, sampleRegion, this::processRecord);

        return mReadCount > 0 && mMaxReadCount > 0;
    }

    private void processRecord(final SAMRecord record)
    {
        ++mReadCount;

        mMaxReadLength = max(mMaxReadLength, record.getReadBases().length);

        if(record.getReadPairedFlag())
        {
            mReadsPaired = true;
            mMateCigarSet |= record.hasAttribute(MATE_CIGAR_ATTRIBUTE);
        }

        if(mConsumer != null)
            mConsumer.accept(record);

        if(mReadCount >= mMaxReadCount)
        {
            mSlicer.haltProcessing();
        }
    }
}
