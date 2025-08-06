package com.hartwig.hmftools.common.bam.testutilities;

import java.util.function.Function;

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
    private int chromosomeIndex;

    public BamRegionWriter(int chromosome, final int readsStart, final int readsStop)
    {
        super("chr" + (chromosome + 1), readsStart, readsStop);
        this.chromosomeIndex = chromosome;
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

    public void writeEntries(SAMFileWriter bamWriter, Function<ChromosomeWindow, Pair<BaseRegion, BaseRegion>> basesFromWindow)
    {
        var window = new ChromosomeWindow(chromosomeIndex, start(), start() + readLength);
        var readNumber = 0;
        while(window.start() < end())
        {
            for(int i = 0; i < depthAtEachStep; i++)
            {
                Pair<BaseRegion, BaseRegion> baseRegionPair = basesFromWindow.apply(window);
                var readName = "A:B:C:" + readNumber;

                Pair<SAMRecord, SAMRecord> reads = new PairedRecordsBuilder(readName, bamWriter.getFileHeader()).build(baseRegionPair);
                bamWriter.addAlignment(reads.getLeft());
                bamWriter.addAlignment(reads.getRight());

                readNumber++;
            }
            window = window.next(stepLength);
        }
    }
}
