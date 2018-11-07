package com.hartwig.hmftools.svanalysis;

import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_STRING;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.List;

import com.google.common.collect.Lists;
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
    private static String COL_CLUSTER_DESC = "ClusterDesc";
    private static String COL_RESOLVED_TYPE = "ResolvedType";
    private static String COL_CLUSTER_COUNT = "ClusterCount";
    private static String COL_CHAIN_ID = "ChainId";

    private List<List<String>> mSourceFileData; // the values to check against
    private List<List<String>> mValidateFileData; // the values to check

    private static String SOURCE_FILENAME = "sv_check_src_filename";
    private static String VALIDATE_FILENAME = "sv_check_val_filename";
    private static String SAMPLE_ID = "sv_check_sample";
    private static String CHECK_ENABLED = "sv_check_enabled";

    private boolean mChecksEnabled;
    private String mSourceFileName;
    private String mValidateFileName;
    private List<String> mSamplesList;

    private BufferedWriter mMismatchWriter;
    private String mOutputDir;

    private List<String> mIdColumns; // by name
    private List<String> mCheckColumns; // by name
    private List<Integer> mSourceIdColumns;
    private List<Integer> mValidateIdColumns;
    private List<Integer> mSourceCheckColumns; // by index into the source file
    private List<Integer> mValidateCheckColumns; // by index into the validation file


    public ResultsChecker()
    {
        mSourceIdColumns = Lists.newArrayList();
        mValidateIdColumns = Lists.newArrayList();
        mIdColumns = Lists.newArrayList();
        mCheckColumns= Lists.newArrayList();
        mSourceCheckColumns = Lists.newArrayList();
        mValidateCheckColumns = Lists.newArrayList();
        mSourceFileData = Lists.newArrayList();
        mValidateFileData = Lists.newArrayList();
        mSamplesList = Lists.newArrayList();
        mLogMismatches = false;
        mAllRecordsMatch = false;
        mBreakOnMismatch = true;
        mChecksEnabled = false;
        mOutputDir = "";
        mMismatchWriter = null;
    }

    public void setLogMismatches(boolean toggle) { mLogMismatches = toggle; }
    public void setBreakOnMismatch(boolean toggle) { mBreakOnMismatch = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(VALIDATE_FILENAME, true, "File to validate");
        options.addOption(SOURCE_FILENAME, true, "File to check against");
        options.addOption(SAMPLE_ID, true, "Optional - specific sample");
        options.addOption(CHECK_ENABLED, true, "Enabled=ON, Disable=OFF or not present");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> samplesList, final String outputDir)
    {
        mSourceFileName = cmd.getOptionValue(SOURCE_FILENAME);
        mValidateFileName = cmd.getOptionValue(VALIDATE_FILENAME);
        mSamplesList.addAll(samplesList);
        mChecksEnabled = cmd.hasOption(CHECK_ENABLED) ? cmd.getOptionValue(CHECK_ENABLED).equals("ON") : false;
        mOutputDir = outputDir;

        return mChecksEnabled;
    }

    public boolean loadData()
    {
        if(mSourceFileName.isEmpty() || mValidateFileName.isEmpty())
            return false;

        if(!loadDataFile(mSourceFileName, mSourceIdColumns, mSourceCheckColumns, mSourceFileData))
        {
            return false;
        }

        if(!loadDataFile(mSourceFileName, mValidateIdColumns, mValidateCheckColumns, mValidateFileData))
        {
            return false;
        }

        if(mSourceFileData.isEmpty() || mValidateFileData.isEmpty())
        {
            mChecksEnabled = false;
            return false;
        }

        return true;
    }

    private int MIN_COLUMN_COUNT = 5;

    private boolean loadDataFile(final String filename, final List<Integer> idColumns, final List<Integer> checkColumns, final List<List<String>> dataList)
    {
        if (filename == null || filename.isEmpty())
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty data CSV file({})", filename);
                return false;
            }

            int invalidRowCount = 0;

            String[] fieldNames = line.split(",");
            int fieldCount = fieldNames.length;

            // check for required columns and extract their indices for use during the matching routine
            List<String> idColumnsNames = Lists.newArrayList();
            idColumnsNames.addAll(mIdColumns);
            List<String> checkColumnsNames = Lists.newArrayList();
            checkColumnsNames.addAll(mCheckColumns);

            for(int i = 0; i < fieldCount; ++i)
            {
                final String fieldName = fieldNames[i];

                if(idColumnsNames.contains(fieldName))
                {
                    idColumnsNames.remove(fieldName);
                    idColumns.add(i);
                }
                else if(checkColumnsNames.contains(fieldName))
                {
                    checkColumnsNames.remove(fieldName);
                    checkColumns.add(i);
                }
            }

            if(!checkColumnsNames.isEmpty() || !idColumnsNames.isEmpty())
            {
                LOGGER.error("not all required columns were matched");
                return false;
            }

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",", -1);

                if(items.length < MIN_COLUMN_COUNT)
                {
                    ++invalidRowCount;
                    continue;
                }

                if(!mSamplesList.isEmpty())
                {
                    final String sampleId = items[COL_SAMPLE_ID_INDEX];
                    if(!mSamplesList.contains(sampleId))
                        continue;
                }

                List<String> dataValues = Lists.newArrayList();

                for (int i = 0; i < items.length; ++i)
                {
                    dataValues.add(items[i]);
                }

                dataList.add(dataValues);
            }

            if(invalidRowCount > 0)
            {
                LOGGER.warn("loaded {} data sets, invalid rows({})", dataList.size(), invalidRowCount);
                return false;
            }

            LOGGER.debug("loaded {} data sets", dataList.size());
            return true;
        }
        catch (IOException exception)
        {
            LOGGER.error("failed to read data file({})", filename);
            return false;
        }
    }

    public void setIdColumns(boolean matchOnPosition)
    {
        mIdColumns.add(COL_SAMPLE_ID);

        if(matchOnPosition)
        {
            mIdColumns.add(COL_CHR_START);
            mIdColumns.add(COL_CHR_END);
            mIdColumns.add(COL_POS_START);
            mIdColumns.add(COL_POS_END);
            mIdColumns.add(COL_ORIENT_START);
            mIdColumns.add(COL_ORIENT_END);
            mIdColumns.add(COL_TYPE);
        }
        else
        {
            mIdColumns.add(COL_VAR_ID);
        }
    }

    public void addDefaultColumnsToCheck()
    {
        mCheckColumns.add(COL_CLUSTER_DESC);
        mCheckColumns.add(COL_RESOLVED_TYPE);
        // columns.add(COL_CLUSTER_ID);
        // columns.add(COL_CLUSTER_COUNT);
        mCheckColumns.add(COL_CHAIN_ID);

        // other candidates:
        // AsmbMatchStart,AsmbMatchEnd
        // FoldbackLnkStart,FoldbackLenStart,FoldbackLnkEnd,FoldbackLenEnd
        // LEStart,LEEnd - for Suspect Line elements
        // LnkSvStart,LnkTypeStart,LnkLenStart,LnkInfoStart,LnkSvEnd,LnkTypeEnd,LnkLenEnd,LnkInfoEnd
    }

    public boolean runChecks()
    {
        // walk through the 2 datasets together, matching each SV on either ID or positional data
        // and the checking that the required fields match, logging any discrepancies
        if(mCheckColumns.isEmpty())
            return true;

        mAllRecordsMatch = true;

        int sourceIndex = 0;
        int validateIndex = 0;
        int mismatchCount = 0;

        while(sourceIndex < mSourceFileData.size())
        {
            if(validateIndex >= mValidateFileData.size())
                break;

            List<String> sourceItems = mSourceFileData.get(sourceIndex);
            List<String> validateItems = null;

            // check that the next validate record matches and if not find it
            int currentValidateIndex = validateIndex;
            boolean matchFound = false;
            boolean indexReset = (validateIndex == 0);
            while(currentValidateIndex < mValidateFileData.size())
            {
                validateItems = mSourceFileData.get(currentValidateIndex);

                if(recordsMatch(sourceItems, validateItems))
                {
                    matchFound = true;
                    break;
                }

                ++currentValidateIndex;

                // break once back around to the start
                if(indexReset && currentValidateIndex == validateIndex)
                    break;

                if(currentValidateIndex >= mValidateFileData.size() && !indexReset)
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

            if(!checkFieldMatches(sourceItems, validateItems))
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
            LOGGER.warn("mismatchCount({}) v {} source records", mismatchCount, mSourceFileData.size());

        try
        {
            if (mMismatchWriter != null)
                mMismatchWriter.close();
        }
        catch(IOException e) { }


        return mAllRecordsMatch;
    }

    private boolean checkFieldMatches(final List<String> sourceItems, final List<String> validateItems)
    {
        for(int i = 0; i < mCheckColumns.size(); ++i)
        {
            int sourceIndex = mSourceCheckColumns.get(i);
            int validateIndex = mValidateCheckColumns.get(i);

            if(!sourceItems.get(sourceIndex).equals(validateItems.get(validateIndex)))
            {
                final String fieldName = mCheckColumns.get(i);

                writeMismatchData(sourceItems, validateItems, fieldName, sourceItems.get(sourceIndex), validateItems.get(validateIndex));

                if(mLogMismatches)
                {
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

    private void writeMismatchData(final List<String> sourceItems, final List<String> validateItems,
            final String fieldName, final String sourceValue, final String validateValue)
    {
        try
        {
            BufferedWriter writer = null;

            if(mMismatchWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mMismatchWriter;
            }
            else
            {
                String outputFileName = mOutputDir;

                if(!outputFileName.endsWith("/"))
                    outputFileName += "/";

                outputFileName += "sva_result_mismatches.csv";

                Path outputFile = Paths.get(outputFileName);

                writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                mMismatchWriter = writer;

                // definitional fields
                writer.write("SampleId,SourceId,CheckId");

                if(mIdColumns.size() > 1)
                {
                    for(final String idCol : mIdColumns)
                        writer.write(String.format(",%s", idCol));
                }

                writer.write(",Field,SourceValue,CheckValue");
                writer.newLine();
            }

            writer.write(String.format(",%s,%s,%s",
                    sourceItems.get(COL_SAMPLE_ID_INDEX), sourceItems.get(COL_VAR_ID_INDEX), validateItems.get(COL_VAR_ID_INDEX)));

            if(mIdColumns.size() > 1)
            {
                for(int idColIndex : mSourceIdColumns)
                    writer.write(String.format(",%s", sourceItems.get(idColIndex)));
            }

            writer.write(String.format(",%s,%s,%s", fieldName, sourceValue, validateValue));
            writer.newLine();

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }

    }

}
