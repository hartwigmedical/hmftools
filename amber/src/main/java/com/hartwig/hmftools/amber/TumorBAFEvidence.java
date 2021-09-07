package com.hartwig.hmftools.amber;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableTumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.collection.Multimaps;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorBAFEvidence implements Callable<TumorBAFEvidence>
{
    private final String mChromosome;
    private final String mBamFile;
    private final TumorBAFFactory mBafFactory;
    private final List<ModifiableTumorBAF> mEvidence;
    private final SamReaderFactory mSamReaderFactory;
    private final GenomePositionSelector<ModifiableTumorBAF> mSelector;
    private final BamSlicer mBamSlicer;
    private final List<GenomeRegion> mBafRegions;

    public TumorBAFEvidence(
            int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig, final String bamFile,
            final SamReaderFactory samReaderFactory, final List<BaseDepth> baseDepths)
    {
        mBafFactory = new TumorBAFFactory(minBaseQuality);
        mChromosome = contig;
        mBamFile = bamFile;
        mSamReaderFactory = samReaderFactory;

        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(typicalReadDepth);
        for(BaseDepth bafRegion : baseDepths)
        {
            builder.addPosition(bafRegion);
        }

        mBafRegions = builder.build();
        mEvidence = baseDepths.stream().map(TumorBAFFactory::create).collect(Collectors.toList());
        mSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(mEvidence));
        mBamSlicer = new BamSlicer(minMappingQuality);
    }

    @NotNull
    public String contig()
    {
        return mChromosome;
    }

    @NotNull
    public List<TumorBAF> evidence()
    {
        return mEvidence.stream().filter(x -> x.tumorIndelCount() == 0).collect(Collectors.toList());
    }

    @Override
    public TumorBAFEvidence call() throws Exception
    {
        try(SamReader reader = mSamReaderFactory.open(new File(mBamFile)))
        {
            mBamSlicer.sliceNoDups(reader, mBafRegions, this::processRecord);
        }

        return this;
    }

    private void processRecord(@NotNull final SAMRecord record)
    {
        mSelector.select(asRegion(record), bafEvidence -> mBafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record)
    {
        return GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
