package prep;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrRecord;

import org.apache.commons.lang3.tuple.Pair;

import feature.FeatureType;
import feature.FeatureValue;

public class ReduxBqrPrep
{
    private final PrepConfig mConfig;
    private final List<ExtendedBqrRecord> mExtendedBqrRecords = new ArrayList<>();

    private final int HI_QUAL_THRESHOLD = 30;

    private static final String KEY_FLD_READ_TYPE = "readType";
    private static final String KEY_FLD_STANDARD_MUTATION = "standardMutation";
    private static final String KEY_FLD_STANDARD_TRINUC_CONTEXT = "standardTrinucContext";
    private static final String KEY_FLD_ORIGINAL_QUAL = "originalQualBin";

    public ReduxBqrPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private static class ExtendedBqrRecord
    {
        private final ConsensusType ReadType;
        private final String StandardMutation;
        private final String StandardTrinucContext;
        private final byte OriginalQuality;
        private final double RecalibratedQuality;
        private final int Count;

        private ExtendedBqrRecord(BqrRecord bqrRecord)
        {
            ReadType = bqrRecord.Key.ReadType;
            OriginalQuality = bqrRecord.Key.Quality;
            RecalibratedQuality = bqrRecord.RecalibratedQuality;
            Count = bqrRecord.Count;

            char refBase = (char) bqrRecord.Key.Ref;
            char altBase = (char) bqrRecord.Key.Alt;

            char standardRefBase = refBase;
            char standardAltBase = altBase;
            if(refBase == 'G' || refBase == 'A')
            {
                standardRefBase = swapDnaBase(refBase);
                standardAltBase = swapDnaBase(altBase);
            }

            String standardMutation = String.format("%s>%s", standardRefBase, standardAltBase);

            String standardTrinucContext = new String(
                    refBase == standardRefBase
                    ? bqrRecord.Key.TrinucleotideContext
                    : reverseComplementBases(bqrRecord.Key.TrinucleotideContext)
            );

            StandardMutation = standardMutation;
            StandardTrinucContext = standardTrinucContext;
        }

        private String getOriginalQualBin()
        {
            if(OriginalQuality < 20)
                return "0-19";

            else if(OriginalQuality < 30)
                return "20-29";

            else
                return "30+";
        }
    }

    private void loadSnvBqrRecords(String sampleId)
    {
        String filePath = BqrFile.generateFilename(mConfig.getReduxDir(sampleId), sampleId);

        List<BqrRecord> bqrRecords = BqrFile.read(filePath);

        List<BqrRecord> bqrRecordsFiltered = bqrRecords.stream().filter(x -> x.Key.Ref != x.Key.Alt).toList();

        List<ExtendedBqrRecord> extendedBqrRecords = bqrRecordsFiltered.stream().map(ExtendedBqrRecord::new).toList();

        // Sort here to control the eventual order in which features are plotted in R
        List<ExtendedBqrRecord> extendedBqrRecordsSorted = extendedBqrRecords.stream()
                .sorted(Comparator.comparing((ExtendedBqrRecord x) -> x.ReadType)
                        .thenComparing(x -> x.StandardMutation)
                        .thenComparing(x -> x.StandardTrinucContext)
                )
                .toList();

        mExtendedBqrRecords.addAll(extendedBqrRecordsSorted);
    }

    private static Map<String, Double> calcMeanChangeInQualPerGroup(Map<String, List<ExtendedBqrRecord>> bqrRecordGroups)
    {
        Map<String, Double> keyResultsMap = new LinkedHashMap<>();

        for(String key : bqrRecordGroups.keySet())
        {
            List<ExtendedBqrRecord> bqrRecordsInGroup = bqrRecordGroups.get(key);

            double totalCountInGroup = bqrRecordsInGroup.stream().mapToDouble(x -> x.Count).sum();

            double weightMeanChangeInQual = 0;
            for(ExtendedBqrRecord bqrRecordInGroup : bqrRecordsInGroup)
            {
                double changeInQual = bqrRecordInGroup.RecalibratedQuality - bqrRecordInGroup.OriginalQuality;
                double weight = bqrRecordInGroup.Count / totalCountInGroup;
                double weightedChangeInQual = changeInQual * weight;

                weightMeanChangeInQual += weightedChangeInQual;
            }

            keyResultsMap.put(key, weightMeanChangeInQual);
        }

        return keyResultsMap;
    }

    private List<FeatureValue<Double>> calcChangeInQualPerTrinucContext()
    {
        List<ExtendedBqrRecord> bqrRecords = mExtendedBqrRecords.stream()
                .filter(x -> x.OriginalQuality >= HI_QUAL_THRESHOLD)
                .toList();

        Map<String, List<ExtendedBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(ExtendedBqrRecord bqrRecord : bqrRecords)
        {
            String key = FeatureValue.keyFromPairs(
                    Pair.of(KEY_FLD_READ_TYPE, bqrRecord.ReadType.toString()),
                    Pair.of(KEY_FLD_STANDARD_MUTATION, bqrRecord.StandardMutation),
                    Pair.of(KEY_FLD_STANDARD_TRINUC_CONTEXT, bqrRecord.StandardTrinucContext)
            );

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        Map<String, Double> meanChangeInQuals = calcMeanChangeInQualPerGroup(bqrRecordGroups);

        return meanChangeInQuals.keySet().stream()
                .map(x -> new FeatureValue<>(x, meanChangeInQuals.get(x), FeatureType.REDUX_BQR_PER_SNV96_CONTEXT))
                .toList();
    }

    private List<FeatureValue<Double>> calcChangeInQualPerOriginalQual()
    {
        Map<String, List<ExtendedBqrRecord>> bqrRecordGroups = new LinkedHashMap<>();
        for(ExtendedBqrRecord bqrRecord : mExtendedBqrRecords)
        {
            String key = FeatureValue.keyFromPairs(
                    Pair.of(KEY_FLD_READ_TYPE, bqrRecord.ReadType.toString()),
                    Pair.of(KEY_FLD_STANDARD_MUTATION, bqrRecord.StandardMutation),
                    Pair.of(KEY_FLD_ORIGINAL_QUAL, bqrRecord.getOriginalQualBin())
            );

            bqrRecordGroups.putIfAbsent(key, new ArrayList<>());
            bqrRecordGroups.get(key).add(bqrRecord);
        }

        Map<String, Double> meanChangeInQuals = calcMeanChangeInQualPerGroup(bqrRecordGroups);

        return meanChangeInQuals.keySet().stream()
                .map(x -> new FeatureValue<>(x, meanChangeInQuals.get(x), FeatureType.REDUX_BQR_PER_ORIG_QUAL))
                .toList();
    }

    public List<FeatureValue<Double>> extractSampleData(String sampleId)
    {
        loadSnvBqrRecords(sampleId);

        List<FeatureValue<Double>> featureValues = new ArrayList<>();

        featureValues.addAll(calcChangeInQualPerOriginalQual());
        featureValues.addAll(calcChangeInQualPerTrinucContext());

        return featureValues;
    }

}
