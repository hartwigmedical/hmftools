package com.hartwig.hmftools.common.basequal.jitter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.sequencing.ConsensusType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// store the data for each MS table for a repeat unit such as (A/T)
// each row is a number of unit repeat in reference
// the row stores the jitter, read counts map.
public class JitterCountsTable
{
    public static final Logger sLogger = LogManager.getLogger(JitterCountsTable.class);

    public class Row
    {
        public final int refNumUnits;
        public int totalReadCount = 0;
        public final Map<Integer, Integer> jitterCounts = new HashMap<>();

        public Row(final int refNumUnits)
        {
            this.refNumUnits = refNumUnits;
        }

        void addRead(int jitter)
        {
            addReads(jitter, 1);
        }

        void addReads(int jitter, int numReads)
        {
            totalReadCount += numReads;
            jitterCounts.merge(jitter, numReads, Integer::sum);
        }

        Set<Integer> getJitterKeys() { return jitterCounts.keySet(); }

        int getJitterReadCount(int jitter)
        {
            return jitterCounts.getOrDefault(jitter, 0);
        }

        void setJitterReadCount(int jitter, int count)
        {
            jitterCounts.put(jitter, count);
        }

        int getTotalReadCount() { return totalReadCount; }

        String getRepeatUnit() { return RepeatUnit; }

        com.hartwig.hmftools.common.sequencing.ConsensusType getConsensusType() { return ConsensusType; }
    }

    public final String RepeatUnit;
    public final ConsensusType ConsensusType;

    // ref num unit to rows
    private final List<Row> mRows = new ArrayList<>();

    JitterCountsTable(final String repeatUnit, final ConsensusType consensusType, final double maxSingleAltSiteContributionPerc)
    {
        RepeatUnit = repeatUnit;
        ConsensusType = consensusType;
    }

    public int repeatUnitLength()
    {
        if(this.RepeatUnit.contains("/"))
            return this.RepeatUnit.split("/")[0].length();
        return 3;
    }

    // summarise the data from
    static JitterCountsTable summariseFrom(
            final String repeatUnit, final ConsensusType consensusType, final double maxSingleAltSiteContributionPerc,
            final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        // In order to filter out outliers, we perform the stats summation in a loop
        // We create a table of all read stats, then use that table to filter out outliers and create a
        // new table. We do this iteratively until no outlier is found.
        JitterCountsTable outlierTestTable = null;
        boolean outlierFound = false;
        Set<MicrosatelliteSiteAnalyser> outliers = Collections.newSetFromMap(new IdentityHashMap<>());

        while(true)
        {
            JitterCountsTable newTable = new JitterCountsTable(repeatUnit, consensusType, maxSingleAltSiteContributionPerc);

            for(MicrosatelliteSiteAnalyser microsatelliteSiteAnalyser : microsatelliteSiteAnalysers)
            {
                if(outliers.contains(microsatelliteSiteAnalyser))
                {
                    continue;
                }

                if(!microsatelliteSiteAnalyser.shouldKeepSite(JitterAnalyserConstants.ALT_COUNT_FRACTION_INIT,
                        JitterAnalyserConstants.ALT_COUNT_FRACTION_STEP,
                        JitterAnalyserConstants.MAX_REJECTED_READ_FRACTION,
                        JitterAnalyserConstants.MIN_PASSING_SITE_READS))
                {
                    continue;
                }

                // get all the read counts into a row object
                Row row = newTable.new Row(microsatelliteSiteAnalyser.refGenomeMicrosatellite().numRepeat);

                for(Map.Entry<Integer, Integer> entry : microsatelliteSiteAnalyser.passingJitterCounts(consensusType).entrySet())
                {
                    int jitter = entry.getKey();
                    int numReads = entry.getValue();
                    row.addReads(jitter, numReads);
                }

                newTable.mergeCounts(row);
            }

            if(outlierTestTable != null && !outlierFound)
            {
                return newTable;
            }

            outlierFound = false;
            outlierTestTable = newTable;
        }
    }

    void mergeCounts(Row rowToMerge)
    {
        Row row = getOrCreateRow(rowToMerge.refNumUnits);

        // add all the reads from this table
        for(Map.Entry<Integer, Integer> entry: rowToMerge.jitterCounts.entrySet())
        {
            row.addReads(entry.getKey(), entry.getValue());
        }
    }

    public List<Row> getRows()
    {
        return mRows;
    }

    public int getReadCount(int refNumRepeats)
    {
        Row row = getRow(refNumRepeats);
        return row == null ? 0 : row.totalReadCount;
    }

    public Row getOrCreateRow(int refNumUnits)
    {
        for(int i = 0; i < mRows.size(); ++i)
        {
            Row row = mRows.get(i);
            if(row.refNumUnits == refNumUnits)
            {
                return row;
            }
            if(row.refNumUnits > refNumUnits)
            {
                // make a row object, insert but keep list sorted
                Row newRow = new Row(refNumUnits);
                mRows.add(i, newRow);
                return newRow;
            }
        }
        // make new object and add to end
        Row row = new Row(refNumUnits);
        mRows.add(row);
        return row;
    }

    @Nullable
    public Row getRow(int refNumUnits)
    {
        for(Row row : mRows)
        {
            if(row.refNumUnits == refNumUnits)
            {
                return row;
            }
        }
        return null;
    }

    public int totalReadCount() { return mRows.stream().mapToInt(Row::getTotalReadCount).sum(); }
}
