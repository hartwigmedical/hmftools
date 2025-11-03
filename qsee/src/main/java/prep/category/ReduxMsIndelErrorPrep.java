package prep.category;

import java.io.File;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterCountsTableFile;
import com.hartwig.hmftools.common.redux.JitterTableRow;

import org.apache.commons.lang3.tuple.Pair;

import feature.Feature;
import feature.FeatureKey;
import feature.FeatureType;
import feature.SourceTool;
import prep.CategoryPrep;
import prep.PrepConfig;

public class ReduxMsIndelErrorPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private static final int MAX_REF_UNIT_COUNT = 15;

    private static final String KEY_FLD_REPEAT_TYPE = "repeatUnitType";
    private static final String KEY_FLD_CONSENSUS_TYPE = "consensusType";
    private static final String KEY_FLD_REF_NUM_UNITS = "refNumUnits";

    public ReduxMsIndelErrorPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private String findBackwardsCompatibleJitterFile(String sampleId) throws NoSuchFileException
    {
        // TODO: Remove this temporary method. In WiGiTS 3.0, the (new) REDUX jitter file path will be used.

        File newJitterFile = new File(JitterCountsTableFile.generateFilename(mConfig.getReduxDir(sampleId), sampleId));
        File oldJitterFile = new File(mConfig.getReduxDir(sampleId) + File.separator + sampleId + ".ms_table.tsv.gz");

        if(newJitterFile.isFile())
        {
            return newJitterFile.getAbsolutePath();
        }

        if(oldJitterFile.isFile())
        {
            return oldJitterFile.getAbsolutePath();
        }

        throw new NoSuchFileException(String.format("sample(%s) - could not determine the REDUX(%s) or SAGE(%s) jitter file",
                sampleId, newJitterFile.getName(), oldJitterFile.getName()));
    }

    private List<JitterTableRow> loadJitterCountsTable(String sampleId) throws IOException
    {
        String filePath = findBackwardsCompatibleJitterFile(sampleId);
        Collection<JitterCountsTable> tables = JitterCountsTableFile.read(filePath);

        List<JitterTableRow> tablesFlattened = new ArrayList<>();
        for(JitterCountsTable table : tables)
        {
            tablesFlattened.addAll(table.getRows());
        }

        List<JitterTableRow> tableFiltered = tablesFlattened.stream()
                .filter(x -> x.refNumUnits() <= MAX_REF_UNIT_COUNT)
                .toList();

        return tableFiltered;
    }

    @VisibleForTesting static final int INSERTIONS_KEY = 1;
    @VisibleForTesting static final int DELETIONS_KEY = -1;

    @VisibleForTesting
    static List<JitterTableRow> aggregateByJitter(List<JitterTableRow> table)
    {
        List<JitterTableRow> newTable = new ArrayList<>();

        for(JitterTableRow row : table)
        {
            int insertionCount = row.jitterCounts().entrySet().stream()
                    .filter(x -> x.getKey() > 0)
                    .mapToInt(x -> x.getValue())
                    .sum();

            int deletionCount = row.jitterCounts().entrySet().stream()
                    .filter(x -> x.getKey() < 0)
                    .mapToInt(x -> x.getValue())
                    .sum();

            JitterTableRow newRow = new JitterTableRow(row.refNumUnits(), row.getRepeatUnit(), row.getConsensusType());
            newRow.setTotalReadCount(row.totalReadCount());
            newRow.setJitterReadCount(INSERTIONS_KEY, insertionCount);
            newRow.setJitterReadCount(DELETIONS_KEY, deletionCount);

            newTable.add(newRow);
        }

        return newTable;
    }

    private static String getRepeatType(String repeatUnit)
    {
        if(repeatUnit.matches("^\\w/\\w$"))
            return repeatUnit + " repeat";

        else if(repeatUnit.matches("^\\w{2}/.*"))
            return "2bp repeat";

        else if(repeatUnit.matches("^\\d+bp repeat"))
            return "≥3bp repeat";

        else
            throw new IllegalArgumentException("Unexpected repeat unit: " + repeatUnit);
    }

    @VisibleForTesting
    static List<JitterTableRow> aggregateByRepeatType(List<JitterTableRow> table)
    {
        Map<FeatureKey, List<JitterTableRow>> groupedTables = new LinkedHashMap<>();
        for(JitterTableRow row : table)
        {
            FeatureKey key = FeatureKey.ofPairs(
                    Pair.of(KEY_FLD_CONSENSUS_TYPE, row.getConsensusType().toString()),
                    Pair.of(KEY_FLD_REPEAT_TYPE, getRepeatType(row.getRepeatUnit())),
                    Pair.of(KEY_FLD_REF_NUM_UNITS, String.valueOf(row.refNumUnits()))
            );

            groupedTables.putIfAbsent(key, new ArrayList<>());
            groupedTables.get(key).add(row);
        }

        List<JitterTableRow> newTable = new ArrayList<>();
        for(FeatureKey key : groupedTables.keySet())
        {
            List<JitterTableRow> groupRows = groupedTables.get(key);

            JitterTableRow oldFirstRow = groupRows.get(0);

            JitterTableRow newRow = new JitterTableRow(
                    oldFirstRow.refNumUnits(),
                    getRepeatType(oldFirstRow.getRepeatUnit()),
                    oldFirstRow.getConsensusType()
            );

            if(groupRows.size() == 1)
            {
                newRow.setTotalReadCount(oldFirstRow.totalReadCount());
                newRow.setJitterReadCount(INSERTIONS_KEY, oldFirstRow.getJitterReadCount(INSERTIONS_KEY));
                newRow.setJitterReadCount(DELETIONS_KEY, oldFirstRow.getJitterReadCount(DELETIONS_KEY));
            }
            else
            {
                int totalReadCount = 0;
                int insertionReadCount = 0;
                int deletionReadCount = 0;
                for(JitterTableRow row : groupRows)
                {
                    totalReadCount += row.totalReadCount();
                    insertionReadCount += row.getJitterReadCount(INSERTIONS_KEY);
                    deletionReadCount += row.getJitterReadCount(DELETIONS_KEY);
                }

                newRow.setTotalReadCount(totalReadCount);
                newRow.setJitterReadCount(INSERTIONS_KEY, insertionReadCount);
                newRow.setJitterReadCount(DELETIONS_KEY, deletionReadCount);
            }

            newTable.add(newRow);
        }

        return newTable;
    }

    private static double calcPhredScore(int count, int totalCount)
    {
        return (totalCount > 0)
                ? -10 * Math.log10(count / (double) totalCount)
                : Double.NaN;
    }

    @VisibleForTesting
    static List<Feature> calcPhredScores(List<JitterTableRow> table)
    {
        List<Feature> indelPhredScores = new ArrayList<>();
        List<Feature> indelPhredScoreDiffs = new ArrayList<>();

        for(JitterTableRow row : table)
        {
            String consensusType = row.getConsensusType().toString();
            String repeatType = row.getRepeatUnit();
            String refNumUnits = String.valueOf(row.refNumUnits());

            int totalReadCount = row.totalReadCount();

            int deletionReadCount = row.getJitterReadCount(DELETIONS_KEY);
            int insertionReadCount = row.getJitterReadCount(INSERTIONS_KEY);
            int indelReadCount = deletionReadCount + insertionReadCount;

            FeatureKey key = FeatureKey.ofPairs(
                    Pair.of(KEY_FLD_CONSENSUS_TYPE, consensusType),
                    Pair.of(KEY_FLD_REPEAT_TYPE, repeatType),
                    Pair.of(KEY_FLD_REF_NUM_UNITS, refNumUnits)
            );

            double indelPhredScore = calcPhredScore(indelReadCount, totalReadCount);
            indelPhredScores.add(new Feature(key.withType(FeatureType.MS_INDEL_ERROR_RATES), indelPhredScore, SourceTool.REDUX));

            double insertionPhredScore = calcPhredScore(insertionReadCount, totalReadCount);
            double deletionPhredScore = calcPhredScore(deletionReadCount, totalReadCount);
            double indelPhredScoreDiff = deletionPhredScore - insertionPhredScore;
            indelPhredScoreDiffs.add(new Feature(key.withType(FeatureType.MS_INDEL_ERROR_BIAS), indelPhredScoreDiff, SourceTool.REDUX));
        }

        List<Feature> features = new ArrayList<>();
        features.addAll(indelPhredScores);
        features.addAll(indelPhredScoreDiffs);

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId) throws IOException
    {
        List<JitterTableRow> table = loadJitterCountsTable(sampleId);

        List<JitterTableRow> tableAggregated;
        tableAggregated = aggregateByJitter(table);
        tableAggregated = aggregateByRepeatType(tableAggregated);

        List<Feature> features = calcPhredScores(tableAggregated);

        return features;
    }
}
