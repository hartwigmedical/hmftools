package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
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
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.prep.bqr.BaseQualBinner;

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

    private List<BqrRecord> loadSnvBqrRecords(String sampleId) throws NoSuchFileException
    {
        String filePath = findBackwardsCompatibleBqrFile(sampleId);

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
    static List<Feature> calcChangeInQualPerTrinucContext(List<BqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        byte hiQualThreshold = baseQualBinner.binRanges().get(BaseQualBin.HIGH).lowerBound();
        bqrRecords = bqrRecords.stream().filter(x -> x.Key.Quality >= hiQualThreshold).toList();

        Map<FeatureKey, List<BqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(BqrRecord bqrRecord : bqrRecords)
        {
            String featureName = FeatureKey.formMultiFieldName(
                    FIELD_READ_TYPE, bqrRecord.Key.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, formMutationString(bqrRecord.Key),
                    FIELD_STANDARD_TRINUC_CONTEXT, new String(bqrRecord.Key.TrinucleotideContext),
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binRangeStringFrom(bqrRecord.Key.Quality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_SNV96_CONTEXT, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    @VisibleForTesting
    static List<Feature> calcChangeInQualPerOriginalQual(List<BqrRecord> bqrRecords, BaseQualBinner baseQualBinner)
    {
        Map<FeatureKey, List<BqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(BqrRecord bqrRecord : bqrRecords)
        {
            String featureName = FeatureKey.formMultiFieldName(
                    FIELD_READ_TYPE, bqrRecord.Key.ReadType.toString(),
                    FIELD_STANDARD_MUTATION, formMutationString(bqrRecord.Key),
                    FIELD_ORIGINAL_QUAL_BIN, baseQualBinner.binRangeStringFrom(bqrRecord.Key.Quality)
            );

            FeatureKey key = new FeatureKey(featureName, FeatureType.BQR_PER_ORIG_QUAL, SOURCE_TOOL);

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        List<Feature> features = calcMeanChangeInQualPerGroup(bqrRecordGroups);
        return features;
    }

    private static List<Feature> calcMeanChangeInQualPerGroup(Map<FeatureKey, List<BqrRecord>> bqrRecordGroups)
    {
        Map<FeatureKey, Double> meanChangeInQuals = new LinkedHashMap<>();
        for(FeatureKey key : bqrRecordGroups.keySet())
        {
            List<BqrRecord> bqrRecordsInGroup = bqrRecordGroups.get(key);

            double totalCountInGroup = bqrRecordsInGroup.stream().mapToDouble(x -> x.Count).sum();

            double weightedMeanedChangeInQual = 0;
            for(BqrRecord bqrRecordInGroup : bqrRecordsInGroup)
            {
                double changeInQual = bqrRecordInGroup.RecalibratedQuality - bqrRecordInGroup.Key.Quality;
                double weight = bqrRecordInGroup.Count / totalCountInGroup;
                double weightedChangeInQual = changeInQual * weight;

                weightedMeanedChangeInQual += weightedChangeInQual;
            }

            meanChangeInQuals.put(key, weightedMeanedChangeInQual);
        }

        List<Feature> features = new ArrayList<>();
        for(FeatureKey key : meanChangeInQuals.keySet())
        {
            features.add(new Feature(key, meanChangeInQuals.get(key)));
        }
        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId) throws NoSuchFileException
    {
        List<BqrRecord> bqrRecords = loadSnvBqrRecords(sampleId);

        List<Feature> features = new ArrayList<>();
        features.addAll(calcChangeInQualPerOriginalQual(bqrRecords, mBaseQualBinner));
        features.addAll(calcChangeInQualPerTrinucContext(bqrRecords, mBaseQualBinner));

        return features;
    }
}