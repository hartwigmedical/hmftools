package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;

import org.junit.Test;

public class GnomadFrequencyAnalysisTest
{
    private final GnomadFrequencySupplier Frequencies = new GnomadFrequencySupplier()
    {
        @Override
        public double getFrequency(String chromosome, int position)
        {
            return position * 1.0 / 1000;
        }
    };

    @Test
    public void getMeanTest()
    {
        List<GenomePosition> data = new ArrayList<>();
        data.add(gp(_1, 100));
        data.add(gp(_2, 100));
        data.add(gp(_2, 400));
        assertEquals(0.2, new GnomadFrequencyAnalysis(Frequencies).getMeanFrequency(data), 0.001);
    }

    private GenomePosition gp(HumanChromosome chromosome, int position)
    {
        return new GenomePositionImpl(chromosome, V38, position);
    }
}
