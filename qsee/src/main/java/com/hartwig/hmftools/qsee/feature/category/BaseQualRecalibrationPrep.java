package com.hartwig.hmftools.qsee.feature.category;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrRecord;

import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.CategoryPrep;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.feature.category.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.feature.category.bqr.BaseQualBinner;
import com.hartwig.hmftools.qsee.feature.category.bqr.StandardBqrRecord;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

public class BaseQualRecalibrationPrep implements CategoryPrep
{
    private final QseeConfig mConfig;
    private final BaseQualBinner mBaseQualBinner;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final String FIELD_READ_TYPE = "ReadType";
    private static final String FIELD_STANDARD_MUTATION = "StandardMutation";
    private static final String FIELD_STANDARD_TRINUC_CONTEXT = "StandardTrinucContext";
    private static final String FIELD_READ_ORIENTATION = "ReadOrientation";
    private static final String FIELD_ORIGINAL_QUAL_BIN = "OriginalQualBin";

    // Field value used when records are not grouped by that field
    private static final String ANY = "ANY";

    /*
    Add a pseudocount to 'count' to avoid divide by zero errors when calculating weighted mean of change in qual.

    The pseudocount is set to a higher value than 1 so that when counts are e.g. [1, 2, 0, 0], the effective weight
    of low counts (e.g. 1 and 2) is reduced. We want to reduce skewed weighted means due to sparse data.
     */
    private static final int PSEUDOCOUNT = 5;

    public BaseQualRecalibrationPrep(QseeConfig config)
    {
        mConfig = config;
        mBaseQualBinner = new BaseQualBinner(config.SequencingTech);
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.BASE_QUAL_RECALIBRATION; }
    
    private List<StandardBqrRecord> loadSnvBqrRecords(String sampleId, SampleType sampleType)
    {
        String baseDir = mConfig.getReduxDir(sampleId, sampleType);
        String filePath = BqrFile.generateFilename(baseDir, sampleId);

        List<BqrRecord> records = BqrFile.read(filePath);
        List<BqrRecord> recordsFiltered = records.stream().filter(x -> x.Key.Ref != x.Key.Alt).toList();
        List<StandardBqrRecord> recordsStandardised = recordsFiltered.stream().map(StandardBqrRecord::from).toList();

        Comparator<StandardBqrRecord> comparator = Comparator.comparing((StandardBqrRecord x) -> x.ReadType())
                .thenComparing(x -> x.StandardContext().mutationString())
                .thenComparing(x -> x.StandardContext().TrinucleotideContext());

        List<StandardBqrRecord> recordsSorted = recordsStandardised.stream().sorted(comparator).toList();

        return recordsSorted;
    }

    private static String formGroupingKey(String readType, String standardMutation, String trinucContext, String readOrientation,
            String originalQualBin)
    {
        return MultiFieldStringBuilder.formMultiField(
                FIELD_READ_TYPE, readType,
                FIELD_STANDARD_MUTATION, standardMutation,
                FIELD_STANDARD_TRINUC_CONTEXT, trinucContext,
                FIELD_READ_ORIENTATION, readOrientation,
                FIELD_ORIGINAL_QUAL_BIN, originalQualBin
        );
    }

    @VisibleForTesting
    static List<Feature> aggregateChangeInQualPerTrinucContext(List<StandardBqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        byte hiQualThreshold = baseQualBinner.binRanges().get(BaseQualBin.HIGH).lowerBound();
        bqrRecords = bqrRecords.stream().filter(x -> x.OriginalQuality() >= hiQualThreshold).toList();

        Map<FeatureKey, List<StandardBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(StandardBqrRecord bqrRecord : bqrRecords)
        {
            String featureName = formGroupingKey(
                    bqrRecord.ReadType().toString(),
                    bqrRecord.StandardContext().mutationString(),
                    bqrRecord.StandardContext().TrinucleotideContext(),
                    bqrRecord.StandardContext().ReadOrientation().name(),
                    baseQualBinner.binRangeStringFrom(bqrRecord.OriginalQuality())
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_SNV96_CONTEXT, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcWeightedMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    @VisibleForTesting
    static List<Feature> aggregateChangeInQualPerOriginalQual(List<StandardBqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        Map<FeatureKey, List<StandardBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(StandardBqrRecord bqrRecord : bqrRecords)
        {
            String featureName = formGroupingKey(
                    bqrRecord.ReadType().toString(),
                    bqrRecord.StandardContext().mutationString(),
                    ANY,
                    bqrRecord.StandardContext().ReadOrientation().name(),
                    baseQualBinner.binRangeStringFrom(bqrRecord.OriginalQuality())
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_ORIG_QUAL, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcWeightedMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    private static List<Feature> calcWeightedMeanChangeInQualPerGroup(Map<FeatureKey, List<StandardBqrRecord>> bqrRecordGroups)
    {
        List<Feature> features = new ArrayList<>();

        for(FeatureKey groupKey : bqrRecordGroups.keySet())
        {
            List<StandardBqrRecord> bqrRecordsInGroup = bqrRecordGroups.get(groupKey);

            double[] changeInQualPerRecord = bqrRecordsInGroup.stream().mapToDouble(x -> x.RecalibratedQuality() - x.OriginalQuality()).toArray();
            double[] variantCountPerRecord = bqrRecordsInGroup.stream().mapToDouble(x -> x.Count() + PSEUDOCOUNT).toArray();

            double weightedMeanChangeInQual = new Mean().evaluate(changeInQualPerRecord, variantCountPerRecord);

            Feature feature = new Feature(groupKey, weightedMeanChangeInQual);
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, SampleType sampleType) throws NoSuchFileException
    {
        List<StandardBqrRecord> bqrRecords = loadSnvBqrRecords(sampleId, sampleType);

        List<Feature> features = new ArrayList<>();
        features.addAll(aggregateChangeInQualPerOriginalQual(bqrRecords, mBaseQualBinner));
        features.addAll(aggregateChangeInQualPerTrinucContext(bqrRecords, mBaseQualBinner));

        return features;
    }
}