package com.hartwig.hmftools.qsee.prep.category;

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

import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class MsIndelErrorPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final int MAX_REF_UNIT_COUNT = 15;

    private static final String FIELD_REPEAT_TYPE = "RepeatUnitType";
    private static final String FIELD_CONSENSUS_TYPE = "ConsensusType";
    private static final String FIELD_REF_NUM_UNITS = "RefNumUnits";

    public MsIndelErrorPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool(){ return SOURCE_TOOL; }

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
        Map<String, List<JitterTableRow>> groupedTables = new LinkedHashMap<>();
        for(JitterTableRow row : table)
        {
            String groupName = MultiFieldStringBuilder.formMultiField(
                    FIELD_CONSENSUS_TYPE, row.getConsensusType().toString(),
                    FIELD_REPEAT_TYPE, getRepeatType(row.getRepeatUnit()),
                    FIELD_REF_NUM_UNITS, String.valueOf(row.refNumUnits())
            );

            groupedTables.putIfAbsent(groupName, new ArrayList<>());
            groupedTables.get(groupName).add(row);
        }

        List<JitterTableRow> newTable = new ArrayList<>();
        for(String groupName : groupedTables.keySet())
        {
            List<JitterTableRow> groupRows = groupedTables.get(groupName);

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
        if(count == 0 || totalCount == 0)
            return Double.NaN;

        return -10 * Math.log10(count / (double) totalCount);
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

            String featureName = MultiFieldStringBuilder.formMultiField(
                    FIELD_CONSENSUS_TYPE, consensusType,
                    FIELD_REPEAT_TYPE, repeatType,
                    FIELD_REF_NUM_UNITS, refNumUnits
            );

            double indelPhredScore = calcPhredScore(indelReadCount, totalReadCount);
            FeatureKey indelPhredKey = new FeatureKey(featureName, FeatureType.MS_INDEL_ERROR_RATES, SOURCE_TOOL);
            indelPhredScores.add(new Feature(indelPhredKey, indelPhredScore));

            double insertionPhredScore = calcPhredScore(insertionReadCount, totalReadCount);
            double deletionPhredScore = calcPhredScore(deletionReadCount, totalReadCount);
            double indelPhredScoreDiff = deletionPhredScore - insertionPhredScore;
            FeatureKey indelPhredScoreDiffKey = new FeatureKey(featureName, FeatureType.MS_INDEL_ERROR_BIAS, SOURCE_TOOL);
            indelPhredScoreDiffs.add(new Feature(indelPhredScoreDiffKey, indelPhredScoreDiff));
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
