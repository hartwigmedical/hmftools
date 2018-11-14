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
    private static final int BACHELOR_CSV_FIELD_COUNT = 17;

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataCollection.class);

    private String mSampleId;
    private List<BachelorGermlineVariant> mGermlineVariants;

    public BachelorDataCollection()
    {
        mSampleId = "";
        mGermlineVariants = Lists.newArrayList();
    }

    public void setSampleId(final String sampleId) {
        mSampleId = sampleId;
    }

    public boolean loadBachelorData(final String filename)
    {
        if (filename.isEmpty())
            return false;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                // CSV fields PATIENT,SOURCE,PROGRAM,ID,GENES,TRANSCRIPT_ID,CHROM,POS,REF,ALTS,EFFECTS

                if (items.length < BACHELOR_CSV_FIELD_COUNT)
                {
                    LOGGER.warn("invalid item count({}), recordIndex({}) in file({})", items.length, mGermlineVariants.size(), filename);
                    return false;
                }

                final String patientId = items[0];

                if (!mSampleId.equals("*") && !mSampleId.equals("") && !mSampleId.contains(patientId)) {
                    continue;
                }

                BachelorGermlineVariant bachRecord = new BachelorGermlineVariant(patientId,
                        items[1],
                        items[2],
                        items[3],
                        items[4],
                        items[5],
                        items[6],
                        Long.parseLong(items[7]),
                        items[8],
                        items[9],
                        items[10],
                        items[11],
                        items[12],
                        Boolean.parseBoolean(items[13]),
                        Integer.parseInt(items[14]),
                        items[15],
                        items[16]);

                mGermlineVariants.add(bachRecord);

                // if(mGermlineVariants.size() > 10000)
                //    break;
            }

            LOGGER.debug("loaded {} bachelor records", mGermlineVariants.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read bachelor input CSV file({})", filename);
            return false;
        }

        return true;
    }

    public final List<BachelorGermlineVariant> getBachelorVariants() {
        return mGermlineVariants;
    }

    private static int FILTER_CSV_FIELD_COUNT = 8;

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

                // CSV fields CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

                if (items.length != FILTER_CSV_FIELD_COUNT)
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
                        items[5],
                        items[6],
                        items[7]);

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

