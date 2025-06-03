package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_BAF;
import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_RATIO;
import static com.hartwig.hmftools.purple.PurpleTestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.PurpleTestUtils.TUMOR_SAMPLE_ID;
import static com.hartwig.hmftools.purple.PurpleTestUtils.buildCobaltChromosomes;
import static com.hartwig.hmftools.purple.PurpleTestUtils.buildDefaultConfigBuilder;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createAmberBaf;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createCobaltRatio;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createObservedRegion;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createSegmentation;
import static com.hartwig.hmftools.purple.FittingConfig.MAX_PLOIDY;
import static com.hartwig.hmftools.purple.FittingConfig.PURITY_INCREMENT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.PurpleTestUtils;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Segmentation;
import com.hartwig.hmftools.purple.somatic.SomaticVariantCache;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

import org.junit.Test;

public class PurityPloidyFitTest
{
    private final PurpleConfig mConfig;
    private final SampleData mSampleData;
    private final AmberData mAmberData;
    private final CobaltData mCobaltData;

    private final SomaticSvCache mSvCache;
    private final SomaticVariantCache mSomaticCache;

    private final RegionFitCalculator mRegionFitCalculator;

    private final ReferenceData mReferenceData;
    private final Segmentation mSegmentation;

    public PurityPloidyFitTest()
    {
        ConfigBuilder configBuilder = buildDefaultConfigBuilder();

        configBuilder.setValue(PURITY_INCREMENT, 0.2);
        configBuilder.setValue(MAX_PLOIDY, 4);
        configBuilder.setValue(PURITY_INCREMENT, 0.2);

        mConfig = PurpleTestUtils.buildPurpleConfig(configBuilder);

        mAmberData = new AmberData(100, Gender.MALE);

        CobaltChromosomes cobaltChromosomes = buildCobaltChromosomes();
        mCobaltData = new CobaltData(cobaltChromosomes);

        mSomaticCache = new SomaticVariantCache(mConfig);
        mSvCache = new SomaticSvCache();

        mRegionFitCalculator = new RegionFitCalculator(cobaltChromosomes, mConfig.Fitting, mAmberData.AverageTumorDepth);

        mReferenceData = new ReferenceData(mConfig);
        mSegmentation = createSegmentation(mReferenceData);

        mSampleData = new SampleData(REF_SAMPLE_ID, TUMOR_SAMPLE_ID, mAmberData, mCobaltData, mSvCache, mSomaticCache);

        buildDefaultProfile();
    }

    private void buildDefaultProfile()
    {
        Chromosome chr1 = HumanChromosome._1;

        // Amber data
        mAmberData.ChromosomeBafs.put(chr1, createAmberBaf(chr1.toString(), 1000, 1.0, 0.5));

        mAmberData.TumorSegments.put(chr1, new PCFPosition(TUMOR_BAF, chr1.toString(), 100000));

        // Cobalt data
        mCobaltData.TumorSegments.put(chr1, new PCFPosition(TUMOR_RATIO, chr1.toString(), 100000));

        List<CobaltRatio> cobaltRatios = Lists.newArrayList();
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 1000, 0.5));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 2000, 0.5));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 3000, 1));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 4000, 0.5));
        mCobaltData.Ratios.put(chr1, cobaltRatios);
    }

    private List<ObservedRegion> buildDefaultObservedRegions()
    {
        List<ObservedRegion> observedRegions = Lists.newArrayList();

        Chromosome chr1 = HumanChromosome._1;

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1, 1000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1001, 2000, 0.5, 0.5, GermlineStatus.DIPLOID, 2));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 2001, 3000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 3001, 4000, 0.5, 0.5, GermlineStatus.DIPLOID, 1));

        return observedRegions;
    }

    @Test
    public void testSomaticFit()
    {
        List<ObservedRegion> observedRegions = buildDefaultObservedRegions();

        PurityPloidyFitter fitter = new PurityPloidyFitter(
                mConfig, mReferenceData, mSampleData, null, mRegionFitCalculator, observedRegions, Gender.FEMALE, false);

        assertTrue(fitter.isValid());

        fitter.run();

        assertNotNull(fitter.finalFit());
    }

    @Test
    public void testMostDiploidPurity()
    {
        FittedPurity fp1 = createFittedPurity(0.3, 0.3, 2.3);
        FittedPurity fp2 = createFittedPurity(0.3, 0.2, 1.9);
        FittedPurity fp3 = createFittedPurity(0.4, 0.4, 1.8);
        FittedPurity fp4 = createFittedPurity(0.4, 0.3, 2.05);

        List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        List<FittedPurity> result = BestFit.mostDiploidPerPurity(all);

        assertEquals(2, result.size());
        assertEquals(fp2, result.get(0));
        assertEquals(fp4, result.get(1));
    }

    private FittedPurity createFittedPurity(double purity, double score, double ploidy)
    {
        return ImmutableFittedPurity.builder()
                .purity(purity)
                .normFactor(0.95)
                .score(score)
                .diploidProportion(1.0)
                .ploidy(ploidy)
                .somaticPenalty(0).build();
    }
}
