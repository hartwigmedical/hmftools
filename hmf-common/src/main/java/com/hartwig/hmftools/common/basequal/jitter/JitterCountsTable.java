package com.hartwig.hmftools.common.basequal.jitter;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.sequencing.ConsensusType;

import org.jetbrains.annotations.Nullable;

public class JitterCountsTable
{
    // store the data for each MS table for a repeat unit such as (A/T)
    // each row is a number of unit repeat in reference
    // the row stores the jitter, read counts map.

    public final String RepeatUnit;
    public final ConsensusType ConsensusType;

    // ref num unit to rows
    private final List<JitterTableRow> mRows = new ArrayList<>();

    public JitterCountsTable(final String repeatUnit, final ConsensusType consensusType, final double maxSingleAltSiteContributionPerc)
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

    public void mergeCounts(JitterTableRow rowToMerge)
    {
        JitterTableRow row = getOrCreateRow(rowToMerge.refNumUnits());

        // add all the reads from this table
        for(Map.Entry<Integer, Integer> entry: rowToMerge.jitterCounts().entrySet())
        {
            row.addReads(entry.getKey(), entry.getValue());
        }
    }

    public List<JitterTableRow> getRows()
    {
        return mRows;
    }

    public int getReadCount(int refNumRepeats)
    {
        JitterTableRow row = getRow(refNumRepeats);
        return row == null ? 0 : row.totalReadCount();
    }

    public JitterTableRow getOrCreateRow(int refNumUnits)
    {
        for(int i = 0; i < mRows.size(); ++i)
        {
            JitterTableRow row = mRows.get(i);
            if(row.refNumUnits() == refNumUnits)
            {
                return row;
            }
            if(row.refNumUnits() > refNumUnits)
            {
                // make a row object, insert but keep list sorted
                JitterTableRow newRow = new JitterTableRow(refNumUnits, RepeatUnit, ConsensusType);
                mRows.add(i, newRow);
                return newRow;
            }
        }
        // make new object and add to end
        JitterTableRow row = new JitterTableRow(refNumUnits, RepeatUnit, ConsensusType);
        mRows.add(row);
        return row;
    }

    @Nullable
    public JitterTableRow getRow(int refNumUnits)
    {
        for(JitterTableRow row : mRows)
        {
            if(row.refNumUnits() == refNumUnits)
            {
                return row;
            }
        }
        return null;
    }

    public int totalReadCount() { return mRows.stream().mapToInt(JitterTableRow::totalReadCount).sum(); }
}
