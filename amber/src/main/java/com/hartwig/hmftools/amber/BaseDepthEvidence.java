package com.hartwig.hmftools.amber;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BaseDepthEvidence implements Callable<BaseDepthEvidence>
{
    private final String mChromosome;
    private final String mBamFile;
    private final SamReaderFactory mSamReaderFactory;
    private final List<ModifiableBaseDepth> EvidenceMap;
    private final GenomePositionSelector<ModifiableBaseDepth> mSelector;
    private final BaseDepthFactory mBafFactory;
    private final SAMSlicer mSamSlicer;

    public BaseDepthEvidence(int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig, final String bamFile,
            final SamReaderFactory samReaderFactory, final List<AmberSite> bafRegions)
    {
        mBafFactory = new BaseDepthFactory(minBaseQuality);
        mChromosome = contig;
        mBamFile = bamFile;
        mSamReaderFactory = samReaderFactory;
        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(typicalReadDepth);
        bafRegions.forEach(builder::addPosition);
        final List<GenomeRegion> bafRegions1 = builder.build();

        EvidenceMap = bafRegions.stream().map(BaseDepthFactory::fromAmberSite).collect(Collectors.toList());
        mSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(EvidenceMap));
        mSamSlicer = new SAMSlicer(minMappingQuality, bafRegions1);
    }

    @NotNull
    public String contig()
    {
        return mChromosome;
    }

    @NotNull
    public List<BaseDepth> evidence()
    {
        return new ArrayList<>(EvidenceMap);
    }

    @Override
    public BaseDepthEvidence call() throws Exception
    {
        try(SamReader reader = mSamReaderFactory.open(new File(mBamFile)))
        {
            mSamSlicer.slice(reader, this::record);
        }

        return this;
    }

    private void record(@NotNull final SAMRecord record)
    {
        mSelector.select(asRegion(record), bafEvidence -> mBafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record)
    {
        return GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
