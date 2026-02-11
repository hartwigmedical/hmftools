package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.*;

import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricSummary;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

import org.junit.Test;

public class SummaryTablePrepTest
{
    @Test
    public void allSummaryTableFeaturesExtracted()
    {
        PurityContext purityContext = createTestPurityContext();
        BamMetricSummary bamMetricSummary = createTestBamMetricSummary();
        BamMetricCoverage bamMetricCoverage = createTestBamMetricCoverage();

        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);
        SummaryTablePrep.putFeatures(purityContext, featuresMap);
        SummaryTablePrep.putFeatures(bamMetricSummary, featuresMap);
        SummaryTablePrep.putFeatures(bamMetricCoverage, featuresMap);

        Set<SummaryTableFeature> allFeatures = EnumSet.allOf(SummaryTableFeature.class);
        Set<SummaryTableFeature> populatedFeatures = featuresMap.keySet();

        assertEquals(SummaryTableFeature.values().length, featuresMap.size());
        assertEquals(allFeatures, populatedFeatures);
    }

    @Test
    public void canCalcPropBasesWithMinCoverage()
    {
        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);
        BamMetricCoverage bamMetricCoverage = createTestBamMetricCoverage();
        SummaryTablePrep.putFeatures(bamMetricCoverage, featuresMap);

        assertEquals(0.9, featuresMap.get(SummaryTableFeature.MIN_COVERAGE_10).value(), 0.01);
        assertEquals(0.7, featuresMap.get(SummaryTableFeature.MIN_COVERAGE_30).value(), 0.01);
        assertEquals(0.6, featuresMap.get(SummaryTableFeature.MIN_COVERAGE_100).value(), 0.01);
        assertEquals(0.4, featuresMap.get(SummaryTableFeature.MIN_COVERAGE_250).value(), 0.01);
    }

    private PurityContext createTestPurityContext()
    {
        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(0.85)
                .ploidy(2.0)
                .normFactor(1.0)
                .score(0.95)
                .diploidProportion(1.0)
                .somaticPenalty(0.0)
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(0.80)
                .maxPurity(0.90)
                .minPloidy(2.0)
                .maxPloidy(2.2)
                .minDiploidProportion(0.75)
                .maxDiploidProportion(0.85)
                .build();

        PurpleQC qc = ImmutablePurpleQC.builder()
                .status(List.of(PurpleQCStatus.PASS))
                .method(FittedPurityMethod.NORMAL)
                .copyNumberSegments(100)
                .unsupportedCopyNumberSegments(1)
                .purity(0.85)
                .amberGender(Gender.FEMALE)
                .cobaltGender(Gender.FEMALE)
                .deletedGenes(10)
                .contamination(0.01)
                .lohPercent(0.01)
                .amberMeanDepth(100)
                .tincLevel(0.1)
                .chimerismPercentage(0.01)
                .build();

        return ImmutablePurityContext.builder()
                .gender(Gender.FEMALE)
                .runMode(RunMode.TUMOR_GERMLINE)
                .targeted(false)
                .bestFit(fittedPurity)
                .method(FittedPurityMethod.NORMAL)
                .score(score)
                .qc(qc)
                .polyClonalProportion(0.1)
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(1.0)
                .tumorMutationalBurdenPerMb(10.0)
                .tumorMutationalLoad(1000)
                .svTumorMutationalBurden(10)
                .microsatelliteStatus(MicrosatelliteStatus.MSS)
                .tumorMutationalLoadStatus(TumorMutationalStatus.HIGH)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.HIGH)
                .build();
    }

    private BamMetricSummary createTestBamMetricSummary()
    {
        return ImmutableBamMetricSummary.builder()
                .totalRegionBases(1_000_000_000L)
                .totalReads(10_000_000L)
                .duplicateReads(1_000_000L)
                .dualStrandReads(1_000_000L)
                .meanCoverage(100.0)
                .sdCoverage(10.0)
                .medianCoverage(100)
                .madCoverage(10)
                .lowMapQualPercent(0.01)
                .duplicatePercent(0.01)
                .unmappedPercent(0.01)
                .lowBaseQualPercent(0.01)
                .overlappingReadPercent(0.01)
                .cappedCoveragePercent(0.01)
                .coverageLevels(List.of(30, 60, 100))
                .coveragePercents(List.of(0.95, 0.90, 0.85))
                .build();
    }

    private BamMetricCoverage createTestBamMetricCoverage()
    {
        List<ValueFrequency> coverageData = List.of(
                new ValueFrequency(5, 1000),
                new ValueFrequency(10, 1000),
                new ValueFrequency(20, 1000),
                new ValueFrequency(50, 1000),
                new ValueFrequency(100, 1000),
                new ValueFrequency(200, 1000),
                new ValueFrequency(500, 1000),
                new ValueFrequency(1000, 1000),
                new ValueFrequency(2000, 1000),
                new ValueFrequency(5000, 1000)
        );

        return new BamMetricCoverage(coverageData);
    }
}
