package com.hartwig.hmftools.bamtools.remapper.testutilities;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamRegionWriter extends ChrBaseRegion
{
    private int readLength = 100;
    private int stepLength = 10;
    private int depthAtEachStep = 1;

    public BamRegionWriter(int chromosome, final int readsStart, final int readsStop)
    {
        super("chr" + (chromosome + 1), readsStart, readsStop);
    }

    public void setReadLength(int readLength)
    {
        this.readLength = readLength;
    }

    public void setStepLength(int stepLength)
    {
        this.stepLength = stepLength;
    }

    public void setDepthAtEachStep(int depthAtEachStep)
    {
        this.depthAtEachStep = depthAtEachStep;
    }

    public void writeEntries(SAMFileWriter bamWriter, RefGenomeSource refGenomeSource)
    {
        var window = new ChromosomeWindow(0, start(), start() + readLength);
        var readNumber = 0;
        while(window.start() < end())
        {
            for(int i = 0; i < depthAtEachStep; i++)
            {
                BaseRegion left = window.toBaseRegion(refGenomeSource);
                BaseRegion right = window.mateBaseRegion(refGenomeSource);
                var readName = "A:B:C:" + readNumber;

                Pair<SAMRecord, SAMRecord> reads = new PairedRecordsBuilder(readName, bamWriter.getFileHeader()).build(left, right);
                bamWriter.addAlignment(reads.getLeft());
                bamWriter.addAlignment(reads.getRight());

                readNumber++;
            }
            window = window.next(stepLength);
        }
    }
}
