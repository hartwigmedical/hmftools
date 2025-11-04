package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.when;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.cobalt.targeted.WholeGenome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.mockito.Mockito;

public class CobaltCalculatorTest extends CalculationsTestBase
{
    // Enrichment factors in targeted mode.
    private static final double EF1 = 1.1;
    private static final double EF2 = 1.2;
    // Tumor read depths, reference read depths,  tumor gc ratios, reference gc ratios for windows on chr 1
    double[] td1 = { -1.0, -1.0, -1.0, 25.0, 25.0, 24.0, 22.0, 23.0, 24.0, 22.0, 26.0, 26.0, 27.0, 23.0, 25.0 };
    double[] rd1 = { -1.0, -1.0, -1.0, 14.0, 12.0, 12.0, 16.0, 16.0, 17.0, 18.0, 18.0, 19.0, 12.0, 12.0, 15.0 };
    double[] tgc1 = { 0.60, 0.40, 0.41, 0.40, 0.40, 0.40, 0.43, 0.42, 0.42, 0.43, 0.41, 0.41, 0.41, 0.42, 0.43 };
    double[] rgc1 = { 0.60, 0.50, 0.51, 0.50, 0.50, 0.50, 0.51, 0.51, 0.51, 0.52, 0.52, 0.52, 0.53, 0.53, 0.53 };
    // Tumor read depths, reference read depths,  tumor gc ratios, reference gc ratios for windows on chr 2
    double[] td2 = { 26.0, 30.0, 30.0, 34.0, 28.0, 32.0 };
    double[] rd2 = { 16.0, 20.0, 18.0, 22.0, 20.0, 24.0 };
    double[] tgc2 = { 0.44, 0.44, 0.46, 0.46, 0.45, 0.45 };
    double[] rgc2 = { 0.54, 0.54, 0.55, 0.55, 0.56, 0.56 };
    Map<Double, Double> refGcRatios = new HashMap<>();
    Map<Double, Double> tumGcRatios = new HashMap<>();
    Map<Double, Double> tumGcRatiosTargeted = new HashMap<>();
    Map<Double, Double> refGcRatiosTargeted = new HashMap<>();
    double TumorMeanMedianRatio;
    double ReferenceMeanMedianRatio;
    CobaltConfig config;

    final ListMultimap<HumanChromosome, DepthReading> tumorDepths = ArrayListMultimap.create();
    final ListMultimap<HumanChromosome, DepthReading> referenceDepths = ArrayListMultimap.create();

    CobaltScope targetedScope = new CobaltScope()
    {
        @Override
        public ResultsNormaliser finalNormaliser()
        {
            return new UnityNormaliser();
        }

        @Override
        public ReadDepthStatisticsNormaliser medianByMeanNormaliser()
        {
            return new NoOpReadDepthStatisticsNormaliser();
        }

        @Override
        public double enrichmentQuotient(final Chromosome chromosome, final int position)
        {
            if(chromosome == _1)
            {
                return EF1;
            }
            return EF2;
        }

        @Override
        public ResultsConsolidator resultsConsolidator(final double medianReadDepth)
        {
            return new NoOpConsolidator();
        }

        @Override
        public boolean onTarget(final Chromosome chromosome, final int position)
        {
            if(chromosome == _1)
            {
                return position > 4000;
            }
            return true;
        }
    };

    public CobaltCalculatorTest()
    {
        // reference values:
        // gc: depths
        // 0.50: 100 (masked out), 14, 12, 12 => median = 12
        // 0.51: 110 (masked out), 16, 16, 17 => median = 16
        // 0.52: 18, 18, 19 => median = 18
        // 0.53: 12, 12, 15 => median = 12
        // 0.54: 16, 20 => median = 18
        // 0.55: 18, 22 => median = 20
        // 0.56: 20, 24 => median = 22
        // Reference gc buckets smoothed average read depths:
        // 49: 12/3, 50: (12+16)/3, 51: (12+16+18)/3, 52: (16+18+12)/3, etc
        refGcRatios.put(0.49, 12.0 / 3.0);
        refGcRatios.put(0.50, 28.0 / 3.0);
        refGcRatios.put(0.51, 46.0 / 3.0);
        refGcRatios.put(0.52, 46.0 / 3.0);
        refGcRatios.put(0.53, 48.0 / 3.0);
        refGcRatios.put(0.54, 50.0 / 3.0);
        refGcRatios.put(0.55, 60.0 / 3.0);
        refGcRatios.put(0.56, 42.0 / 3.0);
        refGcRatios.put(0.57, 22.0 / 3.0);
        // tumour values
        // gc: depths
        // 0.40: 200 (masked out), 25, 25, 24 => median = 25
        // 0.41: 220 (masked out), 26, 26, 27 => median = 26
        // 0.42: 23, 24, 23 => median = 23
        // 0.43: 22, 22, 25 => median = 22
        // 0.44: 26, 30 => median = 28
        // 0.45: 28, 32 => median = 30
        // 0.46: 30, 34 => median = 32
        // Tumor gc buckets smoothed average read depths
        // 39: 25/3, 40: (25+26)/3, 41: (25+26+23)/3, 42: (26+23+22)/3 etc
        tumGcRatios.put(0.39, 25.0 / 3.0);
        tumGcRatios.put(0.40, 51.0 / 3.0);
        tumGcRatios.put(0.41, 74.0 / 3.0);
        tumGcRatios.put(0.42, 71.0 / 3.0);
        tumGcRatios.put(0.43, 73.0 / 3.0);
        tumGcRatios.put(0.44, 80.0 / 3.0);
        tumGcRatios.put(0.45, 90.0 / 3.0);
        tumGcRatios.put(0.46, 62.0 / 3.0);
        tumGcRatios.put(0.47, 32.0 / 3.0);

        // Because the window 1:4001-5000 is not in the targeted regions, it does not contribute
        // to the gc buckets, so these need to be adjusted. The values for 0.40 are now 24 and 25,
        // with a median of 24.5. Only the first few smoothed buckets are affected by this.
        tumGcRatiosTargeted.putAll(tumGcRatios);
        tumGcRatiosTargeted.put(0.39, (24.5) / 3.0);
        tumGcRatiosTargeted.put(0.40, (24.5 + 26) / 3.0);
        tumGcRatiosTargeted.put(0.41, (24.5 + 26 + 23) / 3.0);
        // Similarly for reference values.
        refGcRatiosTargeted.putAll(refGcRatios);
        // 0.50: 100 (masked out), 12, 12 => median = 12
        // 0.51: 110 (masked out), 16, 16, 17 => median = 16
        // 49: 12/3, 50: (12+16)/3, 51: (12+16+18)/3, 52: (16+18+12)/3, etc
        //        refGcRatiosTargeted.put(0.50, (14 + 16) / 3.0);
        //        refGcRatiosTargeted.put(0.51, (24.5 + 26 + 23) / 3.0);

        // Set up gc profile data to have everything mappable but for one region on chr1.
        ListMultimap<Chromosome, GCProfile> gcProfileData = ArrayListMultimap.create();
        int position = 0;
        gcProfileData.put(_1, gcProfile(_1, position, 1.0));
        position += 1000;
        gcProfileData.put(_1, gcProfile(_1, position, 0.80)); // unmappable
        for(int i = 0; i < 13; i++)
        {
            position += 1000;
            gcProfileData.put(_1, gcProfile(_1, position, 1.0));
        }
        position = 0;
        for(int i = 0; i < 6; i++)
        {
            gcProfileData.put(_2, gcProfile(_2, position, 1.0));
            position += 1000;
        }

        // Set up two excluded regions in a single window on chr1
        List<ChrBaseRegion> excludedRegions = new ArrayList<>();
        excludedRegions.add(new ChrBaseRegion(V38.versionedChromosome(_1), 2010, 2160));
        excludedRegions.add(new ChrBaseRegion(V38.versionedChromosome(_1), 2210, 2360));

        // Set up config to return these excluded regions and profiles
        config = Mockito.mock(CobaltConfig.class);
        when(config.gcProfileData()).thenReturn(gcProfileData);
        when(config.excludedRegions()).thenReturn(excludedRegions);
        when(config.diploidRegions()).thenReturn(ArrayListMultimap.create());
        when(config.genome()).thenReturn(V38);

        // Create the reference and tumor depth readings
        for(int i = 0; i < 15; i++)
        {
            position = 1000 * i + 1;
            referenceDepths.put(_1, dr(_1, position, rd1[i], rgc1[i]));
            tumorDepths.put(_1, dr(_1, position, td1[i], tgc1[i]));
        }
        for(int i = 0; i < 6; i++)
        {
            position = 1000 * i + 1;
            referenceDepths.put(_2, dr(_2, position, rd2[i], rgc2[i]));
            tumorDepths.put(_2, dr(_2, position, td2[i], tgc2[i]));
        }

        // Get the median and mean read depths for the included regions.
        DescriptiveStatistics tumorStats = new DescriptiveStatistics();
        DescriptiveStatistics referenceStats = new DescriptiveStatistics();
        for(int i = 3; i < 15; i++)
        {
            tumorStats.addValue(td1[i]);
            referenceStats.addValue(rd1[i]);
        }
        for(int i = 0; i < 6; i++)
        {
            tumorStats.addValue(td2[i]);
            referenceStats.addValue(rd2[i]);
        }
        double tumorMean = tumorStats.getMean();
        double tumorMedian = tumorStats.getPercentile(50);
        TumorMeanMedianRatio = tumorMean / tumorMedian;
        double referenceMean = referenceStats.getMean();
        double referenceMedian = referenceStats.getPercentile(50);
        ReferenceMeanMedianRatio = referenceMean / referenceMedian;
    }

    @Test
    public void tumorAndReferenceWholeGenomeTest()
    {
        when(config.scope()).thenReturn(new WholeGenome());

        CobaltCalculator calculator = new CobaltCalculator(tumorDepths, referenceDepths, config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, rd1[0], td1[0], -1.0, -1.0, -1.0, rgc1[0], tgc1[0]);
        checkRatio(ratios1.get(1), _1, 1001, rd1[1], td1[1], -1.0, -1.0, -1.0, rgc1[1], tgc1[1]);
        checkRatio(ratios1.get(2), _1, 2001, rd1[2], td1[2], -1.0, -1.0, -1.0, rgc1[2], tgc1[2]);
        List<Double> referenceRatios1 = new ArrayList<>();
        for(int i = 3; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double refRatio = rd1[i] / refGcRatios.get(rgc1[i]);
            refRatio = refRatio / ReferenceMeanMedianRatio;
            referenceRatios1.add(refRatio);
            double tumRatio = td1[i] / tumGcRatios.get(tgc1[i]);
            tumRatio = tumRatio / TumorMeanMedianRatio;
            // Diploid normalisation has no effect over these ranges so refRatio == refGcDiploidRatio
            checkRatio(ratios1.get(i), _1, position, rd1[i], td1[i], refRatio, tumRatio, refRatio, rgc1[i], tgc1[i]);
        }

        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        List<Double> referenceRatios2 = new ArrayList<>();
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double refRatio = rd2[i] / refGcRatios.get(rgc2[i]);
            refRatio = refRatio / ReferenceMeanMedianRatio;
            referenceRatios2.add(refRatio);
            double tumRatio = td2[i] / tumGcRatios.get(tgc2[i]);
            tumRatio = tumRatio / TumorMeanMedianRatio;
            checkRatio(ratios2.get(i), _2, position, rd2[i], td2[i], refRatio, tumRatio, refRatio, rgc2[i], tgc2[i]);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(2, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(Doubles.median(referenceRatios1), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals(12, medianRatios.get(0).Count, 0.001);
        assertEquals(V38.versionedChromosome(_2), medianRatios.get(1).Chromosome);
        assertEquals(Doubles.median(referenceRatios2), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals(6, medianRatios.get(1).Count, 0.001);
    }

    @Test
    public void tumorOnlyWholeGenomeTest()
    {
        when(config.scope()).thenReturn(new WholeGenome());

        CobaltCalculator calculator = new CobaltCalculator(tumorDepths, ArrayListMultimap.create(), config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, -1.0, td1[0], -1.0, -1.0, -1.0, -1.0, tgc1[0]);
        checkRatio(ratios1.get(1), _1, 1001, -1.0, td1[1], -1.0, -1.0, -1.0, -1.0, tgc1[1]);
        checkRatio(ratios1.get(2), _1, 2001, -1.0, td1[2], -1.0, -1.0, -1.0, -1.0, tgc1[2]);
        for(int i = 3; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double tumRatio = td1[i] / tumGcRatios.get(tgc1[i]);
            tumRatio = tumRatio / TumorMeanMedianRatio;
            checkRatio(ratios1.get(i), _1, position, -1.0, td1[i], -1.0, tumRatio, -1.0, -1.0, tgc1[i]);
        }

        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double tumRatio = td2[i] / tumGcRatios.get(tgc2[i]);
            tumRatio = tumRatio / TumorMeanMedianRatio;
            checkRatio(ratios2.get(i), _2, position, -1.0, td2[i], -1.0, tumRatio, -1.0, -1.0, tgc2[i]);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(0, medianRatios.size());
    }

    @Test
    public void referenceOnlyWholeGenomeTest()
    {
        when(config.scope()).thenReturn(new WholeGenome());

        CobaltCalculator calculator = new CobaltCalculator(ArrayListMultimap.create(), referenceDepths, config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, rd1[0], -1.0, -1.0, -1.0, -1.0, rgc1[0], -1.0);
        checkRatio(ratios1.get(1), _1, 1001, rd1[1], -1.0, -1.0, -1.0, -1.0, rgc1[1], -1.0);
        checkRatio(ratios1.get(2), _1, 2001, rd1[2], -1.0, -1.0, -1.0, -1.0, rgc1[2], -1.0);
        List<Double> referenceRatios1 = new ArrayList<>();
        for(int i = 3; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double refRatio = rd1[i] / refGcRatios.get(rgc1[i]);
            refRatio = refRatio / ReferenceMeanMedianRatio;
            referenceRatios1.add(refRatio);
            // Diploid normalisation has no effect over these ranges so refRatio == refGcDiploidRatio
            checkRatio(ratios1.get(i), _1, position, rd1[i], -1.0, refRatio, -1.0, refRatio, rgc1[i], -1.0);
        }

        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        List<Double> referenceRatios2 = new ArrayList<>();
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double refRatio = rd2[i] / refGcRatios.get(rgc2[i]);
            refRatio = refRatio / ReferenceMeanMedianRatio;
            referenceRatios2.add(refRatio);
            checkRatio(ratios2.get(i), _2, position, rd2[i], -1.0, refRatio, -1.0, refRatio, rgc2[i], -1.0);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(2, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(Doubles.median(referenceRatios1), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals(12, medianRatios.get(0).Count, 0.001);
        assertEquals(V38.versionedChromosome(_2), medianRatios.get(1).Chromosome);
        assertEquals(Doubles.median(referenceRatios2), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals(6, medianRatios.get(1).Count, 0.001);
    }

    @Test
    public void tumorAndReferenceTargetedTest()
    {
        when(config.scope()).thenReturn(targetedScope);

        CobaltCalculator calculator = new CobaltCalculator(tumorDepths, referenceDepths, config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, rd1[0], td1[0], -1.0, -1.0, -1.0, rgc1[0], tgc1[0]);
        checkRatio(ratios1.get(1), _1, 1001, rd1[1], td1[1], -1.0, -1.0, -1.0, rgc1[1], tgc1[1]);
        checkRatio(ratios1.get(2), _1, 2001, rd1[2], td1[2], -1.0, -1.0, -1.0, rgc1[2], tgc1[2]);
        checkRatio(ratios1.get(3), _1, 3001, rd1[3], td1[3], -1.0, -1.0, -1.0, rgc1[3], tgc1[3]);
        List<Double> referenceRatios1 = new ArrayList<>();
        // In targeted mode the results are normalised so that their mean is 1.0.
        List<Double> expectedRawTumorRatios = new ArrayList<>();
        List<Double> expectedRawReferenceRatios = new ArrayList<>();
        for(int i = 4; i < 15; i++)
        {
            expectedRawTumorRatios.add(td1[i] / (tumGcRatiosTargeted.get(tgc1[i]) * EF1));
            expectedRawReferenceRatios.add(rd1[i] / (refGcRatiosTargeted.get(rgc1[i]) * EF1));
        }
        for(int i = 0; i < 6; i++)
        {
            expectedRawTumorRatios.add(td2[i] / (tumGcRatiosTargeted.get(tgc2[i]) * EF2));
            expectedRawReferenceRatios.add(rd2[i] / (refGcRatiosTargeted.get(rgc2[i]) * EF2));
        }
        double tumorNormalisationFactor = expectedRawTumorRatios.stream().mapToDouble(Double::doubleValue).average().orElseThrow();
        double referenceNormalisationFactor = expectedRawReferenceRatios.stream().mapToDouble(Double::doubleValue).average().orElseThrow();
        for(int i = 4; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double rawRefRatio = rd1[i] / (refGcRatios.get(rgc1[i]) * EF1);
            referenceRatios1.add(rawRefRatio);
            double refRatio = rawRefRatio / referenceNormalisationFactor;
            double tumRatio = td1[i] / (tumGcRatiosTargeted.get(tgc1[i]) * EF1 * tumorNormalisationFactor);
            checkRatio(ratios1.get(i), _1, position, rd1[i], td1[i], refRatio, tumRatio, refRatio, rgc1[i], tgc1[i]);
        }
        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        List<Double> referenceRatios2 = new ArrayList<>();
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double rawRefRatio = rd2[i] / (refGcRatios.get(rgc2[i]) * EF2);
            referenceRatios2.add(rawRefRatio);
            double refRatio = rawRefRatio / referenceNormalisationFactor;
            double tumRatio = td2[i] / (tumGcRatiosTargeted.get(tgc2[i]) * EF2 * tumorNormalisationFactor);
            checkRatio(ratios2.get(i), _2, position, rd2[i], td2[i], refRatio, tumRatio, refRatio, rgc2[i], tgc2[i]);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(2, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(Doubles.median(referenceRatios1), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals(11, medianRatios.get(0).Count, 0.001);
        assertEquals(V38.versionedChromosome(_2), medianRatios.get(1).Chromosome);
        assertEquals(Doubles.median(referenceRatios2), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals(6, medianRatios.get(1).Count, 0.001);
    }

    @Test
    public void tumorOnlyTargetedTest()
    {
        when(config.scope()).thenReturn(targetedScope);

        CobaltCalculator calculator = new CobaltCalculator(tumorDepths, ArrayListMultimap.create(), config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, -1.0, td1[0], -1.0, -1.0, -1.0, -1.0, tgc1[0]);
        checkRatio(ratios1.get(1), _1, 1001, -1.0, td1[1], -1.0, -1.0, -1.0, -1.0, tgc1[1]);
        checkRatio(ratios1.get(2), _1, 2001, -1.0, td1[2], -1.0, -1.0, -1.0, -1.0, tgc1[2]);
        checkRatio(ratios1.get(3), _1, 3001, -1.0, td1[3], -1.0, -1.0, -1.0, -1.0, tgc1[3]);

        // In targeted mode the results are normalised so that their mean is 1.0.
        List<Double> expectedRawTumorRatios = new ArrayList<>();
        for(int i = 4; i < 15; i++)
        {
            expectedRawTumorRatios.add(td1[i] / (tumGcRatiosTargeted.get(tgc1[i]) * EF1));
        }
        for(int i = 0; i < 6; i++)
        {
            expectedRawTumorRatios.add(td2[i] / (tumGcRatiosTargeted.get(tgc2[i]) * EF2));
        }
        double normalisationFactor = expectedRawTumorRatios.stream().mapToDouble(Double::doubleValue).average().orElseThrow();
        for(int i = 4; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double tumRatio = td1[i] / (tumGcRatiosTargeted.get(tgc1[i]) * EF1 * normalisationFactor);
            checkRatio(ratios1.get(i), _1, position, -1.0, td1[i], -1.0, tumRatio, -1.0, -1.0, tgc1[i]);
        }
        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double tumRatio = td2[i] / (tumGcRatiosTargeted.get(tgc2[i]) * EF2 * normalisationFactor);
            checkRatio(ratios2.get(i), _2, position, -1.0, td2[i], -1.0, tumRatio, -1.0, -1.0, tgc2[i]);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(0, medianRatios.size());
    }

    @Test
    public void referenceOnlyTargetedTest()
    {
        when(config.scope()).thenReturn(targetedScope);

        CobaltCalculator calculator = new CobaltCalculator(ArrayListMultimap.create(), referenceDepths, config);
        ListMultimap<Chromosome, CobaltRatio> cobaltRatios = calculator.getCalculatedRatios();
        assertEquals(2, cobaltRatios.keySet().size());
        List<Double> expectedRawReferenceRatios = new ArrayList<>();
        for(int i = 4; i < 15; i++)
        {
            expectedRawReferenceRatios.add(rd1[i] / (refGcRatiosTargeted.get(rgc1[i]) * EF1));
        }
        for(int i = 0; i < 6; i++)
        {
            expectedRawReferenceRatios.add(rd2[i] / (refGcRatiosTargeted.get(rgc2[i]) * EF2));
        }
        double referenceNormalisationFactor = expectedRawReferenceRatios.stream().mapToDouble(Double::doubleValue).average().orElseThrow();
        List<CobaltRatio> ratios1 = cobaltRatios.get(_1);
        assertEquals(15, ratios1.size());
        checkRatio(ratios1.get(0), _1, 1, rd1[0], -1.0, -1.0, -1.0, -1.0, rgc1[0], -1.0);
        checkRatio(ratios1.get(1), _1, 1001, rd1[1], -1.0, -1.0, -1.0, -1.0, rgc1[1], -1.0);
        checkRatio(ratios1.get(2), _1, 2001, rd1[2], -1.0, -1.0, -1.0, -1.0, rgc1[2], -1.0);
        checkRatio(ratios1.get(3), _1, 3001, rd1[3], -1.0, -1.0, -1.0, -1.0, rgc1[3], -1.0);
        List<Double> referenceRatios1 = new ArrayList<>();
        for(int i = 4; i < 15; i++)
        {
            int position = i * 1000 + 1;
            double rawRefRatio = rd1[i] / (refGcRatios.get(rgc1[i]) * EF1);
            referenceRatios1.add(rawRefRatio);
            double refRatio = rawRefRatio / referenceNormalisationFactor;
            referenceRatios1.add(rawRefRatio);
            checkRatio(ratios1.get(i), _1, position, rd1[i], -1.0, refRatio, -1.0, refRatio, rgc1[i], -1.0);
        }
        List<CobaltRatio> ratios2 = cobaltRatios.get(_2);
        assertEquals(6, ratios2.size());
        List<Double> referenceRatios2 = new ArrayList<>();
        for(int i = 0; i < 6; i++)
        {
            int position = i * 1000 + 1;
            double rawRefRatio = rd2[i] / (refGcRatios.get(rgc2[i]) * EF2);
            referenceRatios2.add(rawRefRatio);
            double refRatio = rawRefRatio / referenceNormalisationFactor;
            referenceRatios2.add(rawRefRatio);
            checkRatio(ratios2.get(i), _2, position, rd2[i], -1.0, refRatio, -1.0, refRatio, rgc2[i], -1.0);
        }

        List<MedianRatio> medianRatios = calculator.medianRatios();
        assertEquals(2, medianRatios.size());
        assertEquals(V38.versionedChromosome(_1), medianRatios.get(0).Chromosome);
        assertEquals(Doubles.median(referenceRatios1), medianRatios.get(0).MedianRatio, 0.001);
        assertEquals(11, medianRatios.get(0).Count, 0.001);
        assertEquals(V38.versionedChromosome(_2), medianRatios.get(1).Chromosome);
        assertEquals(Doubles.median(referenceRatios2), medianRatios.get(1).MedianRatio, 0.001);
        assertEquals(6, medianRatios.get(1).Count, 0.001);
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
}
