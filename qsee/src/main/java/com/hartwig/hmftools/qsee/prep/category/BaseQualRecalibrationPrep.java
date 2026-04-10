package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;
import com.hartwig.hmftools.qsee.prep.category.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.prep.category.bqr.BaseQualBinner;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

public class BaseQualRecalibrationPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;
    private final BaseQualBinner mBaseQualBinner;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final String FIELD_READ_TYPE = "ReadType";
    private static final String FIELD_STANDARD_MUTATION = "StandardMutation";
    private static final String FIELD_STANDARD_TRINUC_CONTEXT = "StandardTrinucContext";
    private static final String FIELD_ORIGINAL_QUAL_BIN = "OriginalQualBin";

    /*
    Add a pseudocount to 'count' to avoid divide by zero errors when calculating weighted mean of change in qual.

    The pseudocount is set to a higher value than 1 so that when counts are e.g. [1, 2, 0, 0], the effective weight
    of low counts (e.g. 1 and 2) is reduced. We want to reduce skewed weighted means due to sparse data.
     */
    private static final int PSEUDOCOUNT = 5;

    public BaseQualRecalibrationPrep(QseePrepConfig config)
    {
        mConfig = config;
        mBaseQualBinner = new BaseQualBinner(config.SequencingTech);
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.BASE_QUAL_RECALIBRATION; }
    
    @VisibleForTesting
    static BqrRecord standardiseBases(BqrRecord bqrRecord)
    {
        BqrKey key = bqrRecord.Key;

        byte refBase = key.Ref;
        byte altBase = key.Alt;

        byte standardRefBase = refBase;
        byte standardAltBase = altBase;

        if(refBase == 'G' || refBase == 'A')
        {
            standardRefBase = swapDnaBase(refBase);
            standardAltBase = swapDnaBase(altBase);
        }

        byte[] standardTrinucleotideContext = (refBase == standardRefBase)
                ? key.TrinucleotideContext
                : reverseComplementBases(key.TrinucleotideContext);

        BqrKey standardisedKey = new BqrKey(standardRefBase, standardAltBase, standardTrinucleotideContext,
                bqrRecord.Key.Quality, bqrRecord.Key.ReadType);

        return new BqrRecord(standardisedKey, bqrRecord.Count, bqrRecord.RecalibratedQuality);
    }

    private static String formMutationString(BqrKey key) { return String.format("%c>%c", key.Ref, key.Alt); }

    private List<BqrRecord> loadSnvBqrRecords(String sampleId, SampleType sampleType)
    {
        String baseDir = mConfig.getReduxDir(sampleId, sampleType);
        String filePath = BqrFile.generateFilename(baseDir, sampleId);

        List<BqrRecord> records = BqrFile.read(filePath);
        List<BqrRecord> recordsFiltered = records.stream().filter(x -> x.Key.Ref != x.Key.Alt).toList();
        List<BqrRecord> recordsStandardised = recordsFiltered.stream().map(BaseQualRecalibrationPrep::standardiseBases).toList();

        Comparator<BqrRecord> comparator = Comparator.comparing((BqrRecord x) -> x.Key.ReadType)
                .thenComparing(x -> x.Key.Ref)
                .thenComparing(x -> x.Key.Alt)
                .thenComparing(x -> new String(x.Key.TrinucleotideContext));

        List<BqrRecord> recordsSorted = recordsStandardised.stream().sorted(comparator).toList();

        return recordsSorted;
    }

    @VisibleForTesting
    static List<Feature> aggregateChangeInQualPerTrinucContext(List<BqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        byte hiQualThreshold = baseQualBinner.binRanges().get(BaseQualBin.HIGH).lowerBound();
        bqrRecords = bqrRecords.stream().filter(x -> x.Key.Quality >= hiQualThreshold).toList();

        Map<FeatureKey, List<BqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(BqrRecord bqrRecord : bqrRecords)
        {
            String featureName = MultiFieldStringBuilder.formMultiField(
                    FIELD_READ_TYPE, bqrRecord.Key.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, formMutationString(bqrRecord.Key),
                    FIELD_STANDARD_TRINUC_CONTEXT, new String(bqrRecord.Key.TrinucleotideContext),
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binRangeStringFrom(bqrRecord.Key.Quality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_SNV96_CONTEXT, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcWeightedMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    @VisibleForTesting
    static List<Feature> aggregateChangeInQualPerOriginalQual(List<BqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        Map<FeatureKey, List<BqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(BqrRecord bqrRecord : bqrRecords)
        {
            String featureName = MultiFieldStringBuilder.formMultiField(
                    FIELD_READ_TYPE, bqrRecord.Key.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, formMutationString(bqrRecord.Key),
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binRangeStringFrom(bqrRecord.Key.Quality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_ORIG_QUAL, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcWeightedMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    private static List<Feature> calcWeightedMeanChangeInQualPerGroup(Map<FeatureKey, List<BqrRecord>> bqrRecordGroups)
    {
        List<Feature> features = new ArrayList<>();

        for(FeatureKey groupKey : bqrRecordGroups.keySet())
        {
            List<BqrRecord> bqrRecordsInGroup = bqrRecordGroups.get(groupKey);

            double[] changeInQualPerRecord = bqrRecordsInGroup.stream().mapToDouble(x -> x.RecalibratedQuality - x.Key.Quality).toArray();
            double[] variantCountPerRecord = bqrRecordsInGroup.stream().mapToDouble(x -> x.Count + PSEUDOCOUNT).toArray();

            double weightedMeanChangeInQual = new Mean().evaluate(changeInQualPerRecord, variantCountPerRecord);

            Feature feature = new Feature(groupKey, weightedMeanChangeInQual);
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, SampleType sampleType) throws NoSuchFileException
    {
        List<BqrRecord> bqrRecords = loadSnvBqrRecords(sampleId, sampleType);

        List<Feature> features = new ArrayList<>();
        features.addAll(aggregateChangeInQualPerOriginalQual(bqrRecords, mBaseQualBinner));
        features.addAll(aggregateChangeInQualPerTrinucContext(bqrRecords, mBaseQualBinner));

        return features;
    }
}