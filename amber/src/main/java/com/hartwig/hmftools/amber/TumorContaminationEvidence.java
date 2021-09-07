package com.hartwig.hmftools.amber;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.samtools.BamSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorContaminationEvidence implements Callable<TumorContaminationEvidence>
{
    private final String mContig;
    private final String mBamFile;
    private final SamReaderFactory mSamReaderFactory;
    private final Map<BaseDepth, ModifiableBaseDepth> mEvidenceMap;
    private final GenomePositionSelector<ModifiableBaseDepth> mSelector;
    private final BaseDepthFactory mBafFactory;
    private final BamSlicer mBamSlicer;
    private final List<GenomeRegion> mBafRegions;

    public TumorContaminationEvidence(int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig,
            final String bamFile, final SamReaderFactory samReaderFactory, final List<BaseDepth> baseDepths)
    {
        mBafFactory = new BaseDepthFactory(minBaseQuality);
        mContig = contig;
        mBamFile = bamFile;
        mSamReaderFactory = samReaderFactory;
        mEvidenceMap = Maps.newHashMap();

        final List<ModifiableBaseDepth> tumorRecords = Lists.newArrayList();
        for(BaseDepth baseDepth : baseDepths)
        {
            ModifiableBaseDepth modifiableBaseDepth = BaseDepthFactory.create(baseDepth);
            mEvidenceMap.put(baseDepth, modifiableBaseDepth);
            tumorRecords.add(modifiableBaseDepth);
        }

        mSelector = GenomePositionSelectorFactory.create(tumorRecords);

        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(typicalReadDepth);
        baseDepths.forEach(builder::addPosition);

        mBafRegions = builder.build();
        mBamSlicer = new BamSlicer(minMappingQuality);
    }

    @NotNull
    public String contig()
    {
        return mContig;
    }

    @NotNull
    public List<TumorContamination> evidence()
    {
        final List<TumorContamination> result = Lists.newArrayList();
        for(final Map.Entry<BaseDepth, ModifiableBaseDepth> entry : mEvidenceMap.entrySet())
        {
            final BaseDepth normal = entry.getKey();
            final BaseDepth tumor = entry.getValue();
            if(tumor.altSupport() != 0)
            {
                result.add(ImmutableTumorContamination.builder().from(normal).normal(normal).tumor(tumor).build());
            }
        }

        return result;
    }

    @Override
    public TumorContaminationEvidence call() throws Exception
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
