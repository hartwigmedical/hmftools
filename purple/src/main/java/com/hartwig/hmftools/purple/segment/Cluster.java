package com.hartwig.hmftools.purple.segment;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

public class Cluster implements GenomeRegion
{
    public final String Chromosome;
    public final int Start;
    public int End;
    public final List<PCFPosition> PcfPositions;
    public final List<SVSegment> Variants;

    public Cluster(final String chromosome, final int start, final int end)
    {
        Chromosome = chromosome;
        Start = start;
        End = end;
        PcfPositions = Lists.newArrayList();
        Variants = Lists.newArrayList();
    }

    public List<GenomePosition> ratios()
    {
        return PcfPositions.stream().filter(x -> !x.Source.equals(PCFSource.TUMOR_BAF)).collect(Collectors.toList());
    }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public int start()
    {
        return Start;
    }

    @Override
    public int end()
    {
        return End;
    }
}
