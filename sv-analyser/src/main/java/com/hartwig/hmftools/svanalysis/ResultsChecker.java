package com.hartwig.hmftools.svanalysis;

import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_STRING;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ResultsChecker
{
    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    private boolean mLogMismatches;
    private boolean mBreakOnMismatch;
    private boolean mAllRecordsMatch;

    // id columns
    private static String COL_SAMPLE_ID = "SampleId";
    private static String COL_VAR_ID = "Id";
    private static String COL_CHR_START = "ChrStart";
    private static String COL_CHR_END = "ChrEnd";
    private static String COL_POS_START = "PosStart";
    private static String COL_POS_END = "PosEnd";
    private static String COL_TYPE = "Type";
    private static String COL_ORIENT_START = "OrientStart";
    private static String COL_ORIENT_END = "OrientEnd";
    private static int COL_SAMPLE_ID_INDEX = 0;
    private static int COL_VAR_ID_INDEX = 1;

    // default columns to check by name
    private static String COL_CLUSTER_ID = "ClusterId";
    private static String COL_CLUSTER_COUNT = "ClusterCount";
    private static String COL_CHAIN_ID = "ChainId";

    private GenericDataCollection mSourceFile;
    private GenericDataCollection mValidateFile;

    private static String SOURCE_FILENAME = "sv_check_src_filename";
    private static String VALIDATE_FILENAME = "sv_check_val_filename";
    private static String SAMPLE_ID = "sv_check_sample";

    private String mSourceFileName;
    private String mValidateFileName;
    private String mSpecificSampleId;

    private List<Integer> mIdColumns;
    private List<String> mColumnsToCheck; // by name
    private List<Integer> mSourceColumnsToCheck; // by index into the source file
    private List<Integer> mValidateColumnsToCheck; // by index into the validation file


    public ResultsChecker()
    {
        mIdColumns = Lists.newArrayList();
        mColumnsToCheck = Lists.newArrayList();
        mSourceColumnsToCheck = Lists.newArrayList();
        mValidateColumnsToCheck = Lists.newArrayList();
        mLogMismatches = false;
        mAllRecordsMatch = false;
        mBreakOnMismatch = true;
    }

    public void setLogMismatches(boolean toggle) { mLogMismatches = toggle; }
    public void setBreakOnMismatch(boolean toggle) { mBreakOnMismatch = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(VALIDATE_FILENAME, true, "File to validate");
        options.addOption(SOURCE_FILENAME, true, "File to check against");
        options.addOption(SAMPLE_ID, true, "Optional - specific sample");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        mSourceFileName = cmd.getOptionValue(SOURCE_FILENAME);
        mValidateFileName = cmd.getOptionValue(VALIDATE_FILENAME);
        mSpecificSampleId = cmd.hasOption(SAMPLE_ID) ? cmd.getOptionValue(SAMPLE_ID) : "";

        boolean dataLoaded = loadData();

        return dataLoaded;
    }

    public boolean loadData()
    {
        if(mSourceFileName.isEmpty() || mValidateFileName.isEmpty())
            return false;

        mSourceFile = GenericDataLoader.loadFile(mSourceFileName, GD_TYPE_STRING);
        mValidateFile = GenericDataLoader.loadFile(mValidateFileName, GD_TYPE_STRING);

        return mSourceFile != null && mValidateFile != null;
    }

    public void setIdColumns(boolean matchOnPosition)
    {
        List<String> idColumns = Lists.newArrayList();

        idColumns.add(COL_SAMPLE_ID);

        if(matchOnPosition)
        {
            idColumns.add(COL_CHR_START);
            idColumns.add(COL_CHR_END);
            idColumns.add(COL_POS_START);
            idColumns.add(COL_POS_END);
            idColumns.add(COL_ORIENT_START);
            idColumns.add(COL_ORIENT_END);
            idColumns.add(COL_TYPE);
        }
        else
        {
            idColumns.add(COL_VAR_ID);
        }

        // assumption is that both source and validate file match location of ID columns
        setColumnIndices(mSourceFile, idColumns, mIdColumns);
    }

    public void addDefaultColumnsToCheck()
    {
        List<String> columns = Lists.newArrayList();
        columns.add(COL_CLUSTER_ID);
        columns.add(COL_CLUSTER_COUNT);
        columns.add(COL_CHAIN_ID);

        setColumnsToCheck(columns);
    }

    public boolean runChecks()
    {
        if(mColumnsToCheck.isEmpty())
            return true;

        mAllRecordsMatch = true;

        final List<List<String>> sourceData = mSourceFile.getStringData();
        final List<List<String>> validateData = mValidateFile.getStringData();

        int sourceIndex = 0;
        int validateIndex = 0;

        int mismatchCount = 0;

        while(sourceIndex < sourceData.size())
        {
            if(validateIndex >= validateData.size())
                break;

            List<String> sourceItems = sourceData.get(sourceIndex);
            List<String> validateItems = null;

            // check that the next validate record matches and if not find it
            int currentValidateIndex = validateIndex;
            boolean matchFound = false;
            boolean indexReset = (validateIndex == 0);
            while(currentValidateIndex < validateData.size())
            {
                validateItems = sourceData.get(currentValidateIndex);

                if(recordsMatch(sourceItems, validateItems))
                {
                    matchFound = true;
                    break;
                }

                ++currentValidateIndex;

                // break once back around to the start
                if(indexReset && currentValidateIndex == validateIndex)
                    break;

                if(currentValidateIndex >= validateData.size() && !indexReset)
                {
                    // reset the search index to zero and continuing look for a record match
                    currentValidateIndex = 0;
                    indexReset = true;
                }
            }

            if(!matchFound)
            {
                LOGGER.debug("source record({}) not matched", sourceIndex);

                if(mBreakOnMismatch)
                    break;

                ++sourceIndex;
                continue;
            }

            validateIndex = currentValidateIndex;

            if(!checkColumMatches(sourceItems, validateItems))
            {
                mAllRecordsMatch = false;

                if(mBreakOnMismatch)
                    break;
            }

            ++sourceIndex;
            ++validateIndex;
        }

        if(mAllRecordsMatch)
            LOGGER.debug("all results successfully validated");
        else
            LOGGER.warn("mismatchCount({}) v {} source records", mismatchCount, sourceData.size());

        return mAllRecordsMatch;
    }

    private boolean checkColumMatches(final List<String> sourceItems, final List<String> validateItems)
    {
        for(int i = 0; i < mColumnsToCheck.size(); ++i)
        {
            int sourceIndex = mSourceColumnsToCheck.get(i);
            int validateIndex = mValidateColumnsToCheck.get(i);

            if(!sourceItems.get(sourceIndex).equals(validateItems.get(validateIndex)))
            {
                if(mLogMismatches)
                {
                    final String fieldName = mColumnsToCheck.get(i);

                    LOGGER.debug("mismatch: sample({}) var({}) field({}) values({} vs {})",
                            sourceItems.get(COL_SAMPLE_ID_INDEX), sourceItems.get(COL_VAR_ID_INDEX),
                            fieldName, sourceItems.get(sourceIndex), validateItems.get(validateIndex));
                }
                return false;
            }
        }

        return true;
    }

    private boolean recordsMatch(final List<String> sourceItems, final List<String> validateItems)
    {
        for(int i = 0; i < mIdColumns.size(); ++i)
        {
            if(!sourceItems.get(i).equals(validateItems.get(i)))
                return false;
        }

        return true;
    }

    private void setColumnsToCheck(final List<String> columns)
    {
        mColumnsToCheck.clear();
        mColumnsToCheck.addAll(columns);

        setColumnIndices(mSourceFile, mColumnsToCheck, mSourceColumnsToCheck);
        setColumnIndices(mValidateFile, mColumnsToCheck, mValidateColumnsToCheck);
    }

    private void setColumnIndices(final GenericDataCollection dataCollection, final List<String> columns, List<Integer> columnIndices)
    {
        columnIndices.clear();

        final List<String> fieldNames = dataCollection.getFieldNames();

        for(final String column : columns)
        {
            boolean found = false;
            for(int i = 0; i < fieldNames.size(); ++i)
            {
                if(fieldNames.get(i).equals(column))
                {
                    columnIndices.add(i);
                    found = true;
                    break;
                }
            }

            if(!found)
                columnIndices.add(-1);
        }
    }

}
