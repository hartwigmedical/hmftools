package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.ArrayListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class DiploidRatioLoader implements Consumer<Locatable>
{
    private final ArrayListMultimap<Chromosome, ReadRatio> mResult = ArrayListMultimap.create();
    private final List<ReadRatio> mContigResult = new ArrayList<>();
    private final Collection<Chromosome> mChromosomeList;

    private String mChromosome = null;
    private int mStart = 0;

    public DiploidRatioLoader(final Collection<Chromosome> chromosomes)
    {
        mChromosomeList = chromosomes;
    }

    public DiploidRatioLoader(final Collection<Chromosome> chromosomes, final String diploidBedPath) throws IOException
    {
        this(chromosomes);
        List<BEDFeature> bedFeatures = new ArrayList<>();

        CB_LOGGER.info("Reading diploid regions from {}", diploidBedPath);
        try(final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(diploidBedPath,
                new BEDCodec(),
                false))
        {
            for(BEDFeature bedFeature : reader.iterator())
            {
                bedFeatures.add(bedFeature);
            }
        }

        bedFeatures.forEach(this);
    }

    @Override
    public void accept(@NotNull Locatable bed)
    {
        if(mChromosome == null || !bed.getContig().equals(mChromosome))
        {
            finaliseCurrent();
            mChromosome = bed.getContig();
            mStart = 1;
        }

        createRatio(bed.getContig(), bed.getStart(), bed.getEnd());
        mStart = bed.getEnd() + 1;
    }

    private void createRatio(String contig, int start, int end)
    {
        int position = start;
        while(position < end)
        {
            mContigResult.add(create(contig, position));
            position += WINDOW_SIZE;
        }
    }

    private static ReadRatio create(String contig, int position)
    {
        return ImmutableReadRatio.builder().chromosome(contig).position(position).ratio(1).build();
    }

    private void finaliseCurrent()
    {
        if(mChromosome != null && mStart > 0)
        {
            for (Chromosome c : mChromosomeList)
            {
                if (mChromosome.equals(c.contig))
                {
                    mResult.putAll(c, mContigResult);
                    break;
                }
            }
        }

        mContigResult.clear();
    }

    @NotNull
    public ArrayListMultimap<Chromosome, ReadRatio> build()
    {
        finaliseCurrent();
        return mResult;
    }
}
