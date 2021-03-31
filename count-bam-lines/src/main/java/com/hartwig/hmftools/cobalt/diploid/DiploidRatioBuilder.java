package com.hartwig.hmftools.cobalt.diploid;

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

public class DiploidRatioBuilder implements Consumer<Locatable> {

    private static final int WINDOW_SIZE = 1000;

    private final ListMultimap<Chromosome, ReadRatio> result = ArrayListMultimap.create();
    private final List<ReadRatio> contigResult = Lists.newArrayList();

    private String contig = "";
    private int start = 0;

    public DiploidRatioBuilder() {
    }

    public DiploidRatioBuilder(@NotNull List<BEDFeature> bedFeatures) {
        bedFeatures.forEach(this);
    }

    @Override
    public void accept(@NotNull Locatable bed) {
        if (!bed.getContig().equals(contig)) {
            finaliseCurrent();
            contig = bed.getContig();
            start = 1;
        }

        createRatio(bed.getContig(), bed.getStart(), bed.getEnd());
        start = bed.getEnd() + 1;
    }

    private void createRatio(String contig, int start, int end) {
        int position = start;
        while (position < end) {
            contigResult.add(create(contig, position));
            position += WINDOW_SIZE;
        }
    }

    @NotNull
    private static ReadRatio create(String contig, int position) {
        return ImmutableReadRatio.builder().chromosome(contig).position(position).ratio(1).build();
    }

    private void finaliseCurrent() {
        if (start > 0) {
            result.putAll(HumanChromosome.fromString(contig), contigResult);
        }

        contigResult.clear();
    }

    @NotNull
    public ListMultimap<Chromosome, ReadRatio> build() {
        finaliseCurrent();
        return result;
    }
}
