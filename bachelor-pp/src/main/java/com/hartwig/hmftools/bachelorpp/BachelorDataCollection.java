package com.hartwig.hmftools.bachelorpp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelorpp.types.BachelorRecordFilter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BachelorDataCollection
{
    private static final int BACHELOR_CSV_FIELD_COUNT = 16;

    private static final int COL_INDEX_SAMPLE = 0;
    private static final int COL_INDEX_SOURCE = 1;
    private static final int COL_INDEX_PROGRAM = 2;
    private static final int COL_INDEX_SV_ID = 3;
    private static final int COL_INDEX_GENE = 4;
    private static final int COL_INDEX_TRAN_ID = 5;
    private static final int COL_INDEX_CHR = 6;
    private static final int COL_INDEX_POS = 7;
    private static final int COL_INDEX_REF = 8;
    private static final int COL_INDEX_ALT = 9;
    private static final int COL_INDEX_EFFECTS = 10;
    private static final int COL_INDEX_ANNOTS = 11;
    private static final int COL_INDEX_PROTEIN = 12;
    private static final int COL_INDEX_HZ = 13;
    private static final int COL_INDEX_PHRED = 14;
    private static final int COL_INDEX_CODING = 15;
    private static final int COL_INDEX_MATCH_TYPE = 16;
    private static final int COL_INDEX_GL_ALT_COUNT = 17;
    private static final int COL_INDEX_GL_READ_DEPTH = 18;
    private static final int COL_INDEX_TUMOR_ALT_COUNT = 19;
    private static final int COL_INDEX_TUMOR_READ_DEPTH = 20;
    private static final int COL_INDEX_CODON_INFO = 21;

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataCollection.class);

    private String mSampleId;
    private List<String> mLimitedSampleList;
    private List<BachelorGermlineVariant> mGermlineVariants;
    private int mFileIndex;
    private int mMaxReadCount;

    private BufferedReader mFileReader;

    public BachelorDataCollection()
    {
        mSampleId = "";
        mGermlineVariants = Lists.newArrayList();
        mFileIndex = 0;
        mMaxReadCount = 0;
        mFileReader = null;
    }

    public void setSampleId(final String sampleId, final List<String> limitedSampleList)
    {
        mSampleId = sampleId;
        mLimitedSampleList = limitedSampleList;
    }
    public void setMaxReadCount(int maxReadCount) { mMaxReadCount = maxReadCount; }

    public boolean loadBachelorData(final String filename)
    {
        if (filename.isEmpty())
            return false;

        try
        {
            mFileReader = new BufferedReader(new FileReader(filename));

            mFileReader.readLine(); // skip header
        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read bachelor input CSV file({})", filename);
            return false;
        }

        return true;
    }

    public boolean processBachelorData()
    {
        mGermlineVariants.clear();

        try
        {
            String line = null;
            boolean exitNextSample = false;
            String currentSampleId = "";

            while ((line = mFileReader.readLine()) != null)
            {
                if(line.isEmpty())
                    break;

                ++mFileIndex;

                // parse CSV data
                String[] items = line.split(",");

                if (items.length < BACHELOR_CSV_FIELD_COUNT)
                {
                    LOGGER.warn("sample({}) invalid item count({}), fileIndex({})",
                            currentSampleId, items.length, mFileIndex);
                    continue;
                }

                final String sampleId = items[COL_INDEX_SAMPLE];

                if (!mSampleId.equals("*") && !mSampleId.equals("") && !mSampleId.contains(sampleId))
                    continue;

                if(!mLimitedSampleList.isEmpty() && !mLimitedSampleList.contains(sampleId))
                    continue;

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;

                    if(exitNextSample)
                    {
                        LOGGER.info("halting read at {} records, fileIndex({})", mGermlineVariants.size(), mFileIndex);
                        return true;
                    }
                }

                // check for annotations with ',' which impact string splitting
                if(items.length > COL_INDEX_CODON_INFO+1)
                {
                    checkAnnotationItems(items);

                    // other error checking
                    boolean duplicateFieldFound = false;
                    for(int j = COL_INDEX_SOURCE+1; j < items.length; ++j)
                    {
                        if(items[j].equals(items[COL_INDEX_SOURCE]))
                        {
                            duplicateFieldFound = true;
                            break;
                        }
                    }

                    if(duplicateFieldFound)
                    {
                        LOGGER.warn("sample({}) skipping duplicated record at fileIndex({})", sampleId, mFileIndex);
                        continue;
                    }
                }

                try
                {

                    // extra fields from newer versions
                    final String matchType = items.length > COL_INDEX_MATCH_TYPE ? items[COL_INDEX_MATCH_TYPE] : "";
                    final String codonInfo = items.length > COL_INDEX_CODON_INFO ? items[COL_INDEX_CODON_INFO] : "";

                    BachelorGermlineVariant bachRecord = new BachelorGermlineVariant(sampleId,
                            items[COL_INDEX_SOURCE],
                            items[COL_INDEX_PROGRAM],
                            items[COL_INDEX_SV_ID],
                            items[COL_INDEX_GENE],
                            items[COL_INDEX_TRAN_ID],
                            items[COL_INDEX_CHR],
                            Long.parseLong(items[COL_INDEX_POS]),
                            items[COL_INDEX_REF],
                            items[COL_INDEX_ALT],
                            items[COL_INDEX_EFFECTS],
                            items[COL_INDEX_ANNOTS],
                            items[COL_INDEX_PROTEIN],
                            Boolean.parseBoolean(items[COL_INDEX_HZ]),
                            Integer.parseInt(items[COL_INDEX_PHRED]),
                            items[COL_INDEX_CODING],
                            matchType,
                            codonInfo);

                    if (items.length > COL_INDEX_TUMOR_READ_DEPTH)
                    {
                        int glAltCount = Integer.parseInt(items[COL_INDEX_GL_ALT_COUNT]);
                        int glReadDepth = Integer.parseInt(items[COL_INDEX_GL_READ_DEPTH]);
                        int tumorAltCount = Integer.parseInt(items[COL_INDEX_TUMOR_ALT_COUNT]);
                        int tumorReadDepth = Integer.parseInt(items[COL_INDEX_TUMOR_READ_DEPTH]);

                        bachRecord.setReadData(glAltCount, glReadDepth, tumorAltCount, tumorReadDepth);
                    }

                    mGermlineVariants.add(bachRecord);
                }
                catch(Exception nfe)
                {
                    LOGGER.debug("line parse error({}) fileIndex({}) line: {}", nfe.toString(), mFileIndex, line);
                    return false;
                }

                if(!exitNextSample && mMaxReadCount > 0 && mGermlineVariants.size() >= mMaxReadCount)
                    exitNextSample = true;
            }

            LOGGER.debug("loaded {} bachelor records", mGermlineVariants.size());
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read bachelor input CSV file: {}", e.toString());
            return false;
        }

        return !mGermlineVariants.isEmpty();
    }

    private static void checkAnnotationItems(String[] items)
    {
        String extraAnnots = "";
        boolean hasExtraAnnots = false;
        int extraItems = 0;

        for (int i = COL_INDEX_ANNOTS + 1; i < items.length; ++i)
        {
            extraAnnots += "," + items[i];
            ++extraItems;

            if (items[i].contains("||"))
            {
                hasExtraAnnots = true;
                break;
            }
        }

        if (!hasExtraAnnots)
            return;

        // otherwise shift them back down
        items[COL_INDEX_ANNOTS] += extraAnnots;

        for (int i = COL_INDEX_ANNOTS + 1; i <= COL_INDEX_CODON_INFO; ++i)
        {
            items[i] = items[i + extraItems];
        }
    }

    public final List<BachelorGermlineVariant> getBachelorVariants() { return mGermlineVariants; }

    private static int FILTER_CSV_FIELD_COUNT = 26;

    public static List<BachelorRecordFilter> loadBachelorFilters(final String filename)
    {
        if (filename.isEmpty())
            return Lists.newArrayList();

        List<BachelorRecordFilter> filterRecords = Lists.newArrayList();

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split("\t");

                //              0   1   2  3   4   5    6           13         17         25
                // CSV fields CHROM,POS,ID,REF,ALT,QUAL,FILTER, ... CLNDN, .. CLNSIG, ... MC

//                 [1] "CHROM"        "POS"          "ID"           "REF"          "ALT"          "QUAL"         "FILTER"       "AF_ESP"       "AF_EXAC"      "AF_TGP"
//                        [11] "ALLELEID"     "CLNDISDB"     "CLNDISDBINCL" "CLNDN"        "CLNDNINCL"    "CLNHGVS"      "CLNREVSTAT"   "CLNSIG"       "CLNSIGCONF"   "CLNSIGINCL"
//                        [21] "CLNVC"        "CLNVCSO"      "CLNVI"        "DBVARID"      "GENEINFO"     "MC"           "ORIGIN"       "RS"           "SSR"

                if (items.length < FILTER_CSV_FIELD_COUNT)
                {
                    LOGGER.warn("invalid item count({}), recordIndex({}) in file({})", items.length, filterRecords.size(), filename);
                    return filterRecords;
                }

                BachelorRecordFilter filterRecord = new BachelorRecordFilter(
                        items[0],
                        Long.parseLong(items[1]),
                        items[2],
                        items[3],
                        items[4],
                        items[13],
                        items[17],
                        items[25]);

                filterRecords.add(filterRecord);
            }

            LOGGER.debug("loaded {} filter records", filterRecords.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read bachelor input CSV file({})", filename);
            return filterRecords;
        }

        return filterRecords;
    }

}

