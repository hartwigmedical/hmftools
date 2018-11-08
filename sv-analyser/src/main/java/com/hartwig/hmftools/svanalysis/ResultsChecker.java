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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private static int COL_SAMPLE_ID_INDEX = 0;
    private static int COL_VAR_ID_INDEX = 1;
    private static int COL_TYPE_INDEX = 2;
    private static int COL_CHR_START_INDEX = 3;
    private static int COL_POS_START_INDEX = 4;
    private static int COL_ORIENT_START_INDEX = 5;
    private static int COL_CHR_END_INDEX = 6;
    private static int COL_POS_END_INDEX = 7;
    private static int COL_ORIENT_END_INDEX = 8;

    // default columns to check by name
    private static String COL_CLUSTER_ID = "ClusterId";
    private static String COL_CLUSTER_DESC = "ClusterDesc";
    private static String COL_RESOLVED_TYPE = "ResolvedType";
    private static String COL_CLUSTER_COUNT = "ClusterCount";
    private static String COL_CHAIN_ID = "ChainId";

    private Map<String, List<List<String>>> mSampleSourceData; // the values to check against, keyed by sampleId
    private Map<String, List<List<String>>> mSampleValidateData; // the values to check

    private static String SOURCE_FILENAME = "sv_check_src_filename";
    private static String VALIDATE_FILENAME = "sv_check_val_filename";
    private static String SAMPLE_ID = "sv_check_sample";
    private static String CHECK_ENABLED = "sv_check_enabled";
    private static String BREAK_ON_MISMATCH = "sv_check_break_on_mismatch";

    private boolean mChecksEnabled;
    private String mSourceFileName;
    private String mValidateFileName;
    private List<String> mSamplesList;
    private boolean mMatchOnSvIds;

    private BufferedWriter mMismatchWriter;
    private String mOutputDir;

    private List<Integer> mIdColumns; // by name
    private List<String> mCheckColumns; // by name
    private List<Integer> mSourceCheckColumns; // by index into the source file
    private List<Integer> mValidateCheckColumns; // by index into the validation file


    public ResultsChecker()
    {
        mIdColumns = Lists.newArrayList();
        mCheckColumns= Lists.newArrayList();
        mSourceCheckColumns = Lists.newArrayList();
        mValidateCheckColumns = Lists.newArrayList();
        mSampleSourceData = new HashMap();
        mSampleValidateData = new HashMap();
        mSamplesList = Lists.newArrayList();
        mLogMismatches = false;
        mAllRecordsMatch = false;
        mBreakOnMismatch = false;
        mChecksEnabled = false;
        mMatchOnSvIds = false;
        mOutputDir = "";
        mMismatchWriter = null;
    }

    public void setLogMismatches(boolean toggle) { mLogMismatches = toggle; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(VALIDATE_FILENAME, true, "File to validate");
        options.addOption(SOURCE_FILENAME, true, "File to check against");
        options.addOption(SAMPLE_ID, true, "Optional - specific sample");
        options.addOption(CHECK_ENABLED, true, "Enabled=ON, Disable=OFF or not present");
        options.addOption(BREAK_ON_MISMATCH, false, "Exit run on first mismatch");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> samplesList, final String outputDir)
    {
        mSourceFileName = cmd.getOptionValue(SOURCE_FILENAME);

        if(!cmd.hasOption(VALIDATE_FILENAME))
        {
            if(samplesList.size() != 1)
            {
                return false;
            }

            mValidateFileName = outputDir + samplesList.get(0) + ".csv";
        }
        else
        {
            mValidateFileName = cmd.getOptionValue(VALIDATE_FILENAME);
        }

        mSamplesList.addAll(samplesList);
        mChecksEnabled = cmd.hasOption(CHECK_ENABLED) ? cmd.getOptionValue(CHECK_ENABLED).equals("ON") : false;
        mBreakOnMismatch = cmd.hasOption(BREAK_ON_MISMATCH);
        mOutputDir = outputDir;

        if(!samplesList.isEmpty())
        {
            LOGGER.debug("checking {} specific samples", samplesList.size());
        }

        return mChecksEnabled;
    }

    public boolean loadData()
    {
        if(mSourceFileName.isEmpty() || mValidateFileName.isEmpty())
            return false;

        if(!loadDataFile(mSourceFileName, mSourceCheckColumns, mSampleSourceData))
        {
            return false;
        }

        if(!loadDataFile(mValidateFileName, mValidateCheckColumns, mSampleValidateData))
        {
            return false;
        }

        if(mSampleSourceData.isEmpty() || mSampleValidateData.isEmpty())
        {
            mChecksEnabled = false;
            return false;
        }

        return true;
    }

    private int MIN_COLUMN_COUNT = 5;

    private boolean loadDataFile(final String filename, final List<Integer> checkColumns, final Map<String, List<List<String>>> sampleDataMap)
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
            List<String> checkColumnsNames = Lists.newArrayList();
            checkColumnsNames.addAll(mCheckColumns);

            for(int i = 0; i < fieldCount; ++i)
            {
                final String fieldName = fieldNames[i];

                if(checkColumnsNames.contains(fieldName))
                {
                    checkColumnsNames.remove(fieldName);
                    checkColumns.add(i);
                }
            }

            if(!checkColumnsNames.isEmpty())
            {
                LOGGER.error("not all required columns to chekc were found");
                return false;
            }

            String currentSample = "";
            List<List<String>> sampleDataList = null;
            int rowCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",", -1);

                if(items.length < MIN_COLUMN_COUNT)
                {
                    ++invalidRowCount;
                    continue;
                }

                final String sampleId = items[COL_SAMPLE_ID_INDEX];

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    if (!mSamplesList.isEmpty())
                    {
                        if (!mSamplesList.contains(sampleId))
                            continue;
                    }

                    sampleDataList = Lists.newArrayList();
                    sampleDataMap.put(sampleId, sampleDataList);
                    currentSample = sampleId;
                }

                List<String> dataValues = Lists.newArrayList();

                for (int i = 0; i < items.length; ++i)
                {
                    dataValues.add(items[i]);
                }

                sampleDataList.add(dataValues);
                ++rowCount;
            }

            if(invalidRowCount > 0)
            {
                LOGGER.info("loaded {} samples, {} records, invalid rows({})", sampleDataMap.size(), rowCount, invalidRowCount);
                return false;
            }

            LOGGER.warn("loaded {} samples, {} records", sampleDataMap.size(), rowCount);
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
        mIdColumns.add(COL_SAMPLE_ID_INDEX);

        if(matchOnPosition)
        {
            mMatchOnSvIds = false;
            mIdColumns.add(COL_CHR_START_INDEX);
            mIdColumns.add(COL_CHR_END_INDEX);
            mIdColumns.add(COL_POS_START_INDEX);
            mIdColumns.add(COL_POS_END_INDEX);
            mIdColumns.add(COL_ORIENT_START_INDEX);
            mIdColumns.add(COL_ORIENT_END_INDEX);
            mIdColumns.add(COL_TYPE_INDEX);
        }
        else
        {
            mMatchOnSvIds = true;
            mIdColumns.add(COL_VAR_ID_INDEX);
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

        for(final Map.Entry<String, List<List<String>>> entry : mSampleSourceData.entrySet())
        {
            final String sampleId = entry.getKey();
            final List<List<String>> sourceData = entry.getValue();

            final List<List<String>> validateData = mSampleValidateData.get(sampleId);

            if(validateData == null)
            {
                LOGGER.error("sampleId({}) not found in validate data", sampleId);
                mAllRecordsMatch = false;
                continue;
            }

            while (sourceIndex < sourceData.size())
            {
                if (validateIndex >= validateData.size())
                    break;

                List<String> sourceItems = sourceData.get(sourceIndex);
                List<String> validateItems = null;

                // check that the next validate record matches and if not find it
                int currentValidateIndex = validateIndex;
                boolean matchFound = false;
                boolean indexReset = (validateIndex == 0);
                while (currentValidateIndex < validateData.size())
                {
                    validateItems = validateData.get(currentValidateIndex);

                    if (recordsMatch(sourceItems, validateItems))
                    {
                        matchFound = true;
                        break;
                    }

                    ++currentValidateIndex;

                    // break once back around to the start
                    if (indexReset && currentValidateIndex == validateIndex)
                        break;

                    if (currentValidateIndex >= validateData.size() && !indexReset)
                    {
                        // reset the search index to zero and continuing look for a record match
                        currentValidateIndex = 0;
                        indexReset = true;
                    }
                }

                if (!matchFound)
                {
                    LOGGER.debug("source record({}) not matched", sourceIndex);

                    if (mBreakOnMismatch)
                        break;

                    ++sourceIndex;
                    continue;
                }

                validateIndex = currentValidateIndex;

                if (!checkFieldMatches(sourceItems, validateItems))
                {
                    mAllRecordsMatch = false;
                    ++mismatchCount;

                    if (mBreakOnMismatch)
                        break;
                }

                ++sourceIndex;
                ++validateIndex;
            }
        }

        if(mAllRecordsMatch)
            LOGGER.debug("all results successfully validated");
        else
            LOGGER.warn("mismatchCount({})", mismatchCount);

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
                    LOGGER.debug("mismatch: sample({}) var({}) field({}) source({}) vs check({})",
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
                if(mMatchOnSvIds)
                {
                    writer.write("SampleId,SvId");
                }
                else
                {
                    writer.write("SampleId,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,SourceId,CheckId");
                }

                writer.write(",Field,SourceValue,CheckValue");
                writer.newLine();
            }
            writer.write(String.format("%s", sourceItems.get(COL_SAMPLE_ID_INDEX)));

            if(mMatchOnSvIds)
            {
                writer.write(String.format(",%s", sourceItems.get(COL_VAR_ID_INDEX)));
            }
            else
            {
                for(int i = 1; i < mIdColumns.size(); ++i)
                    writer.write(String.format(",%s", sourceItems.get(i)));

                writer.write(String.format(",%s,%s",
                        sourceItems.get(COL_VAR_ID_INDEX), validateItems.get(COL_VAR_ID_INDEX)));
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
