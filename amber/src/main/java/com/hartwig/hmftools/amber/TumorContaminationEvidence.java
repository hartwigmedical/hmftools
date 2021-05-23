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
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorContaminationEvidence implements Callable<TumorContaminationEvidence>
{
    private final String contig;
    private final String bamFile;
    private final SamReaderFactory samReaderFactory;
    private final Map<BaseDepth, ModifiableBaseDepth> evidenceMap;
    private final GenomePositionSelector<ModifiableBaseDepth> selector;
    private final BaseDepthFactory bafFactory;
    private final SAMSlicer supplier;

    public TumorContaminationEvidence(int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig,
            final String bamFile, final SamReaderFactory samReaderFactory, final List<BaseDepth> baseDepths)
    {
        this.bafFactory = new BaseDepthFactory(minBaseQuality);
        this.contig = contig;
        this.bamFile = bamFile;
        this.samReaderFactory = samReaderFactory;
        this.evidenceMap = Maps.newHashMap();

        final List<ModifiableBaseDepth> tumorRecords = Lists.newArrayList();
        for(BaseDepth baseDepth : baseDepths)
        {
            ModifiableBaseDepth modifiableBaseDepth = BaseDepthFactory.create(baseDepth);
            evidenceMap.put(baseDepth, modifiableBaseDepth);
            tumorRecords.add(modifiableBaseDepth);
        }
        this.selector = GenomePositionSelectorFactory.create(tumorRecords);

        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(typicalReadDepth);
        baseDepths.forEach(builder::addPosition);
        this.supplier = new SAMSlicer(minMappingQuality, builder.build());
    }

    @NotNull
    public String contig()
    {
        return contig;
    }

    @NotNull
    public List<TumorContamination> evidence()
    {
        final List<TumorContamination> result = Lists.newArrayList();
        for(final Map.Entry<BaseDepth, ModifiableBaseDepth> entry : evidenceMap.entrySet())
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
        try(SamReader reader = samReaderFactory.open(new File(bamFile)))
        {
            supplier.slice(reader, this::record);
        }

        return this;
    }

    private void record(@NotNull final SAMRecord record)
    {
        selector.select(asRegion(record), bafEvidence -> bafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record)
    {
        return GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
