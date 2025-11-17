package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.sage.SageCommon.SAGE_FILE_ID;

import java.io.File;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrRecord;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.prep.bqr.BaseQualBinner;
import com.hartwig.hmftools.qsee.prep.bqr.ExtendedBqrRecord;

public class BaseQualRecalibrationPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;
    private final BaseQualBinner mBaseQualBinner;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final String FIELD_READ_TYPE = "ReadType";
    private static final String FIELD_STANDARD_MUTATION = "StandardMutation";
    private static final String FIELD_STANDARD_TRINUC_CONTEXT = "StandardTrinucContext";
    private static final String FIELD_ORIGINAL_QUAL_BIN = "OriginalQualBin";

    public BaseQualRecalibrationPrep(CommonPrepConfig config)
    {
        mConfig = config;
        mBaseQualBinner = new BaseQualBinner(config.SEQUENCING_TYPE);
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private String findBackwardsCompatibleBqrFile(String sampleId) throws NoSuchFileException
    {
        // TODO: Remove this temporary method. In WiGiTS 3.0, the (new) REDUX BQR file path will be used.

        File reduxBqrFile = new File(BqrFile.generateFilename(mConfig.getReduxDir(sampleId), sampleId));
        File sageBqrFile = new File(mConfig.getSageDir(sampleId) + File.separator + sampleId + SAGE_FILE_ID + ".bqr.tsv");

        if(reduxBqrFile.isFile())
            return reduxBqrFile.getAbsolutePath();

        if(sageBqrFile.isFile())
            return sageBqrFile.getAbsolutePath();

        throw new NoSuchFileException(reduxBqrFile.getName() + " or " + sageBqrFile.getName());
    }

    private List<ExtendedBqrRecord> loadSnvBqrRecords(String sampleId) throws NoSuchFileException
    {
        String filePath = findBackwardsCompatibleBqrFile(sampleId);

        List<BqrRecord> bqrRecords = BqrFile.read(filePath);
        List<BqrRecord> bqrRecordsFiltered = bqrRecords.stream().filter(x -> x.Key.Ref != x.Key.Alt).toList();
        List<ExtendedBqrRecord> extendedBqrRecords = bqrRecordsFiltered.stream().map(ExtendedBqrRecord::from).toList();

        // Sort here to control the eventual order in which features are plotted in R
        Comparator<ExtendedBqrRecord> comparator = Comparator.comparing((ExtendedBqrRecord x) -> x.ReadType)
                .thenComparing(x -> x.StandardMutation)
                .thenComparing(x -> x.StandardTrinucContext);

        List<ExtendedBqrRecord> extendedBqrRecordsSorted = extendedBqrRecords.stream().sorted(comparator).toList();

        return extendedBqrRecordsSorted;
    }

    private static Map<FeatureKey, Double> calcMeanChangeInQualPerGroup(Map<FeatureKey, List<ExtendedBqrRecord>> bqrRecordGroups)
    {
        Map<FeatureKey, Double> keyResultsMap = new LinkedHashMap<>();

        for(FeatureKey key : bqrRecordGroups.keySet())
        {
            List<ExtendedBqrRecord> bqrRecordsInGroup = bqrRecordGroups.get(key);

            double totalCountInGroup = bqrRecordsInGroup.stream().mapToDouble(x -> x.Count).sum();

            double weightedMeanedChangeInQual = 0;
            for(ExtendedBqrRecord bqrRecordInGroup : bqrRecordsInGroup)
            {
                double changeInQual = bqrRecordInGroup.RecalibratedQuality - bqrRecordInGroup.OriginalQuality;
                double weight = bqrRecordInGroup.Count / totalCountInGroup;
                double weightedChangeInQual = changeInQual * weight;

                weightedMeanedChangeInQual += weightedChangeInQual;
            }

            keyResultsMap.put(key, weightedMeanedChangeInQual);
        }

        return keyResultsMap;
    }

    @VisibleForTesting
    static List<Feature> calcChangeInQualPerTrinucContext(List<ExtendedBqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        byte hiQualThreshold = baseQualBinner.binRanges().get(BaseQualBin.HIGH).lowerBound();

        bqrRecords = bqrRecords.stream()
                .filter(x -> x.OriginalQuality >= hiQualThreshold)
                .toList();

        Map<FeatureKey, List<ExtendedBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(ExtendedBqrRecord bqrRecord : bqrRecords)
        {
            String featureName = FeatureKey.formMultiFieldName(
                    FIELD_READ_TYPE, bqrRecord.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, bqrRecord.StandardMutation,
                    FIELD_STANDARD_TRINUC_CONTEXT, bqrRecord.StandardTrinucContext,
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binNameFrom(bqrRecord.OriginalQuality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_SNV96_CONTEXT, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        Map<FeatureKey, Double> meanChangeInQuals = calcMeanChangeInQualPerGroup(bqrRecordGroups);

        return meanChangeInQuals.keySet().stream()
                .map(x -> new Feature(x, meanChangeInQuals.get(x)))
                .toList();
    }

    @VisibleForTesting
    static List<Feature> calcChangeInQualPerOriginalQual(List<ExtendedBqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        Map<FeatureKey, List<ExtendedBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(ExtendedBqrRecord bqrRecord : bqrRecords)
        {
            String featureName = FeatureKey.formMultiFieldName(
                    FIELD_READ_TYPE, bqrRecord.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, bqrRecord.StandardMutation,
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binNameFrom(bqrRecord.OriginalQuality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_ORIG_QUAL, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        Map<FeatureKey, Double> meanChangeInQuals = calcMeanChangeInQualPerGroup(bqrRecordGroups);

        return meanChangeInQuals.keySet().stream()
                .map(x -> new Feature(x, meanChangeInQuals.get(x)))
                .toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId) throws NoSuchFileException
    {
        List<ExtendedBqrRecord> extendedBqrRecords = loadSnvBqrRecords(sampleId);

        List<Feature> features = new ArrayList<>();
        features.addAll(calcChangeInQualPerOriginalQual(extendedBqrRecords, mBaseQualBinner));
        features.addAll(calcChangeInQualPerTrinucContext(extendedBqrRecords, mBaseQualBinner));

        return features;
    }
}