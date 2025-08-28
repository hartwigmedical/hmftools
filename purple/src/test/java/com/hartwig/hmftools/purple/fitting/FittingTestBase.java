package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_BAF;
import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_RATIO;
import static com.hartwig.hmftools.purple.FittingConfig.MAX_PLOIDY;
import static com.hartwig.hmftools.purple.FittingConfig.MIN_PURITY;
import static com.hartwig.hmftools.purple.FittingConfig.PURITY_INCREMENT;
import static com.hartwig.hmftools.purple.FittingTestUtils.buildCobaltChromosomes;
import static com.hartwig.hmftools.purple.FittingTestUtils.createAmberBaf;
import static com.hartwig.hmftools.purple.MiscTestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.MiscTestUtils.SAMPLE_ID;
import static com.hartwig.hmftools.purple.MiscTestUtils.buildDefaultConfigBuilder;
import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurpleConfig;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.FittingTestUtils;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.somatic.SomaticVariantCache;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

public class FittingTestBase
{
    PurpleConfig mConfig;
    SampleData mSampleData;
    AmberData mAmberData;
    CobaltData mCobaltData;

    SomaticSvCache mSvCache;
    SomaticVariantCache mSomaticCache;

    RegionFitCalculator mRegionFitCalculator;

    ReferenceData mReferenceData;

    public void setup()
    {
        ConfigBuilder configBuilder = buildDefaultConfigBuilder();

        configBuilder.setValue(MIN_PURITY, 0.2);
        configBuilder.setValue(PURITY_INCREMENT, 0.2);
        configBuilder.setValue(MAX_PLOIDY, 4);
        configBuilder.setValue(PURITY_INCREMENT, 0.2);

        mConfig = buildPurpleConfig(configBuilder);

        mAmberData = new AmberData(100, Gender.MALE);

        CobaltChromosomes cobaltChromosomes = buildCobaltChromosomes();
        mCobaltData = new CobaltData(cobaltChromosomes);

        mSomaticCache = new SomaticVariantCache(mConfig);
        mSvCache = new SomaticSvCache();

        mRegionFitCalculator = new RegionFitCalculator(cobaltChromosomes, mConfig.Fitting, mAmberData.AverageTumorDepth);

        mReferenceData = new ReferenceData();

        mSampleData = new SampleData(REF_SAMPLE_ID, SAMPLE_ID, mAmberData, mCobaltData, mSvCache, mSomaticCache);

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
        cobaltRatios.add(FittingTestUtils.createCobaltRatio(chr1, 1000, 0.5, 0.5));
        cobaltRatios.add(FittingTestUtils.createCobaltRatio(chr1, 2000, 0.5, 0.5));
        cobaltRatios.add(FittingTestUtils.createCobaltRatio(chr1, 3000, 1, 0.5));
        cobaltRatios.add(FittingTestUtils.createCobaltRatio(chr1, 4000, 0.5, 0.5));
        mCobaltData.Ratios.put(chr1, cobaltRatios);
    }

    protected FittedPurity createFittedPurity(double purity, double score, double ploidy)
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
