package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.bed.BEDFeature;

public class DiploidRatioBuilder implements Consumer<Locatable>
{
    private final ListMultimap<Chromosome, ReadRatio> mResult;
    private final List<ReadRatio> mContigResult;

    private String mChromosome;
    private int mStart;

    public DiploidRatioBuilder()
    {
        mResult = ArrayListMultimap.create();
        mContigResult = Lists.newArrayList();
        mChromosome = "";
        mStart = 0;
    }

    public DiploidRatioBuilder(final List<BEDFeature> bedFeatures)
    {
        this();
        bedFeatures.forEach(this);
    }

    @Override
    public void accept(@NotNull Locatable bed)
    {
        if(!bed.getContig().equals(mChromosome))
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
        if(mStart > 0)
        {
            mResult.putAll(HumanChromosome.fromString(mChromosome), mContigResult);
        }

        mContigResult.clear();
    }

    @NotNull
    public ListMultimap<Chromosome, ReadRatio> build()
    {
        finaliseCurrent();
        return mResult;
    }
}
