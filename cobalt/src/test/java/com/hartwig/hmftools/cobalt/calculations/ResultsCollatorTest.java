package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Test;

public class ResultsCollatorTest extends CalculationsTestBase
{
    final ListMultimap<Chromosome, BamRatio> tumorResults = ArrayListMultimap.create();
    final ListMultimap<Chromosome, BamRatio> referenceResults = ArrayListMultimap.create();

    public ResultsCollatorTest()
    {
        Chromosome chromosome = _1;
        int positon = 1;
        tumorResults.put(chromosome, br(chromosome, positon, -1.0, -1.0));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, -1.0));
        positon = 1001;
        tumorResults.put(chromosome, br(chromosome, positon, -1.0, -1.0));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, -1.0));
        positon = 2001;
        tumorResults.put(chromosome, br(chromosome, positon, 20, 0.40));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, -1.0));
        positon = 3001;
        tumorResults.put(chromosome, br(chromosome, positon, 25, 0.40));
        referenceResults.put(chromosome, br(chromosome, positon, 10, 0.40, 0.12));
        positon = 4001;
        tumorResults.put(chromosome, br(chromosome, positon, 26, 0.40));
        referenceResults.put(chromosome, br(chromosome, positon, 9, 0.40, 0.13));
        positon = 5001;
        tumorResults.put(chromosome, br(chromosome, positon, 24, 0.39));
        referenceResults.put(chromosome, br(chromosome, positon, 11, 0.40, 0.10));
        positon = 6001;
        tumorResults.put(chromosome, br(chromosome, positon, -1.0, 0.39));
        referenceResults.put(chromosome, br(chromosome, positon, 15, 0.48, 0.11));

        chromosome = _2;
        positon = 1;
        tumorResults.put(chromosome, br(chromosome, positon, -1.0, -1.0));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, -1.0));
        positon = 1001;
        tumorResults.put(chromosome, br(chromosome, positon, 35, 0.53));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, -1.0));
        positon = 2001;
        tumorResults.put(chromosome, br(chromosome, positon, 32, 0.54));
        referenceResults.put(chromosome, br(chromosome, positon, 21, 0.32, -1.0));
        positon = 3001;
        tumorResults.put(chromosome, br(chromosome, positon, 25, 0.40));
        referenceResults.put(chromosome, br(chromosome, positon, -1.0, -1.0, 0.12));
        positon = 4001;
        tumorResults.put(chromosome, br(chromosome, positon, -1.0, 0.40));
        referenceResults.put(chromosome, br(chromosome, positon, 28, 0.36, 0.28));
    }

    @Test
    public void collateTest()
    {
        ResultsCollator collator = new ResultsCollator(V38);
        ListMultimap<Chromosome, CobaltRatio> collated = collator.collateResults(tumorResults, referenceResults);
        assertEquals(2, collated.asMap().size());
        List<CobaltRatio> ratios1 = collated.get(_1);
        assertEquals(7, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(1), _1, 1001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(2), _1, 2001, -1.0, 20.0, -1.0, 0.20, -1.0, -1.0, 0.40);
        checkRatio(ratios1.get(3), _1, 3001, 10.0, 25.0, 0.10, 0.25, 0.12, 0.40, 0.40);
        checkRatio(ratios1.get(4), _1, 4001, 9.0, 26.0, 0.09, 0.26, 0.13, 0.40, 0.40);
        checkRatio(ratios1.get(5), _1, 5001, 11.0, 24.0, 0.11, 0.24, 0.10, 0.40, 0.39);
        checkRatio(ratios1.get(6), _1, 6001, 15.0, -1.0, 0.15, -1.0, 0.11, 0.48, 0.39);
        List<CobaltRatio> ratios2 = collated.get(_2);
        assertEquals(5, ratios2.size());
        checkRatio(ratios2.get(0), _2, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios2.get(1), _2, 1001, -1.0, 35.0, -1.0, 0.35, -1.0, -1.0, 0.53);
        checkRatio(ratios2.get(2), _2, 2001, 21.0, 32.0, 0.21, 0.32, -1.0, 0.32, 0.54);
        checkRatio(ratios2.get(3), _2, 3001, -1.0, 25.0, -1.0, 0.25, 0.12, -1.0, 0.40);
        checkRatio(ratios2.get(4), _2, 4001, 28.0, -1.0, 0.28, -1.0, 0.28, 0.36, 0.40);
    }

    @Test
    public void referenceOnlyTest()
    {
        ResultsCollator collator = new ResultsCollator(V38);
        ListMultimap<Chromosome, CobaltRatio> collated = collator.collateResults(ArrayListMultimap.create(), referenceResults);
        assertEquals(2, collated.asMap().size());
        List<CobaltRatio> ratios1 = collated.get(_1);
        assertEquals(7, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(1), _1, 1001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(2), _1, 2001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(3), _1, 3001, 10.0, -1.0, 0.10, -1.0, 0.12, 0.40, -1.0);
        checkRatio(ratios1.get(4), _1, 4001, 9.0, -1.0, 0.09, -1.0, 0.13, 0.40, -1.0);
        checkRatio(ratios1.get(5), _1, 5001, 11.0, -1.0, 0.11, -1.0, 0.10, 0.40, -1.0);
        checkRatio(ratios1.get(6), _1, 6001, 15.0, -1.00, 0.15, -1.00, 0.11, 0.48, -1.0);
        List<CobaltRatio> ratios2 = collated.get(_2);
        assertEquals(5, ratios2.size());
        checkRatio(ratios2.get(0), _2, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios2.get(1), _2, 1001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios2.get(2), _2, 2001, 21.0, -1.0, 0.21, -1.0, -1.0, 0.32, -1.0);
        checkRatio(ratios2.get(3), _2, 3001, -1.0, -1.0, -1.0, -1.0, 0.12, -1.0, -1.0);
        checkRatio(ratios2.get(4), _2, 4001, 28.0, -1.00, 0.28, -1.00, 0.28, 0.36, -1.0);
    }

    @Test
    public void tumorOnlyTest()
    {
        ResultsCollator collator = new ResultsCollator(V38);
        ListMultimap<Chromosome, CobaltRatio> collated = collator.collateResults(tumorResults, ArrayListMultimap.create());
        assertEquals(2, collated.asMap().size());
        List<CobaltRatio> ratios1 = collated.get(_1);
        assertEquals(7, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(1), _1, 1001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios1.get(2), _1, 2001, -1.0, 20.0, -1.0, 0.20, -1.0, -1.0, 0.40);
        checkRatio(ratios1.get(3), _1, 3001, -1.0, 25.0, -1.0, 0.25, -1.0, -1.0, 0.40);
        checkRatio(ratios1.get(4), _1, 4001, -1.0, 26.0, -1.0, 0.26, -1.0, -1.0, 0.40);
        checkRatio(ratios1.get(5), _1, 5001, -1.0, 24.0, -1.0, 0.24, -1.0, -1.0, 0.39);
        checkRatio(ratios1.get(6), _1, 6001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.39);
        List<CobaltRatio> ratios2 = collated.get(_2);
        assertEquals(5, ratios2.size());
        checkRatio(ratios2.get(0), _2, 1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
        checkRatio(ratios2.get(1), _2, 1001, -1.00, 35.0, -1.0, 0.35, -1.0, -1.0, 0.53);
        checkRatio(ratios2.get(2), _2, 2001, -1.0, 32.0, -1.0, 0.32, -1.0, -1.0, 0.54);
        checkRatio(ratios2.get(3), _2, 3001, -1.00, 25.0, -1.0, 0.25, -1.0, -1.0, 0.40);
        checkRatio(ratios2.get(4), _2, 4001, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.40);
    }

    private void checkRatio(CobaltRatio ratio, Chromosome chromosome, int position,
            double referenceReadDepth,
            double tumorReadDepth,
            double referenceGcRatio,
            double tumorGcRatio,
            double referenceGcDiploidRatio,
            double referenceGcContent,
            double tumorGcContent)
    {
        assertEquals(V38.versionedChromosome(chromosome), ratio.chromosome());
        assertEquals(position, ratio.position());
        assertEquals(referenceReadDepth, ratio.referenceReadDepth(), 0.0001);
        assertEquals(tumorReadDepth, ratio.tumorReadDepth(), 0.0001);
        assertEquals(referenceGcRatio, ratio.referenceGCRatio(), 0.0001);
        assertEquals(tumorGcRatio, ratio.tumorGCRatio(), 0.0001);
        assertEquals(referenceGcDiploidRatio, ratio.referenceGCDiploidRatio(), 0.0001);
        assertEquals(referenceGcContent, ratio.referenceGcContent(), 0.0001);
        assertEquals(tumorGcContent, ratio.tumorGcContent(), 0.0001);
    }

    private BamRatio br(Chromosome chromosome, int pos, double depth, double gc)
    {
        BamRatio result = br(chromosome, pos, depth, gc, true);
        if(depth > 0)
        {
            result.normaliseByMean(100.0);
        }
        return result;
    }

    private BamRatio br(Chromosome chromosome, int pos, double depth, double gc, double diploidAdjustedRatio)
    {
        BamRatio result = br(chromosome, pos, depth, gc);
        result.setDiploidAdjustedRatio(diploidAdjustedRatio);
        return result;
    }
}
