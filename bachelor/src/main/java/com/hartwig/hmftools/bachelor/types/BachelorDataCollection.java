package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.BACHELOR_CSV_FIELD_COUNT;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.COL_INDEX_ANNOTS;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.COL_INDEX_CLINVAR_SIG_INFO;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.COL_INDEX_SAMPLE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BachelorDataCollection
{
    private static final Logger LOGGER = LogManager.getLogger(BachelorDataCollection.class);

    private List<BachelorGermlineVariant> mGermlineVariants;
    private int mFileIndex;

    private BufferedReader mFileReader;

    public BachelorDataCollection()
    {
        mGermlineVariants = Lists.newArrayList();
        mFileIndex = 0;
        mFileReader = null;
    }

    public boolean loadBachelorData(final String filename)
    {
        mGermlineVariants.clear();

        if (filename.isEmpty())
            return false;

        try
        {
            mFileReader = new BufferedReader(new FileReader(filename));

            mFileReader.readLine(); // skip header

            String line = null;

            while ((line = mFileReader.readLine()) != null)
            {
                if(line.isEmpty())
                    break;

                ++mFileIndex;

                // parse CSV data
                String[] items = line.split(",", -1);

                if (items.length < BACHELOR_CSV_FIELD_COUNT)
                {
                    LOGGER.warn("invalid item count({}), fileIndex({})", items.length, mFileIndex);
                    continue;
                }

                final String sampleId = items[COL_INDEX_SAMPLE];

                // check for annotations with ',' which impacts string splitting
                if(items.length > BACHELOR_CSV_FIELD_COUNT)
                {
                    checkAnnotationItems(items);
                }

                try
                {

                    BachelorGermlineVariant bachRecord = BachelorGermlineVariant.fromCsv(items);
                    mGermlineVariants.add(bachRecord);
                }
                catch(Exception nfe)
                {
                    LOGGER.debug("line parse error({}) fileIndex({}) line: {}", nfe.toString(), mFileIndex, line);
                    return false;
                }
            }

            LOGGER.debug("loaded {} bachelor records", mGermlineVariants.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("failed to read bachelor input file({}): {}", filename, exception.toString());
            return false;
        }

        return true;
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

        for (int i = COL_INDEX_ANNOTS + 1; i <= COL_INDEX_CLINVAR_SIG_INFO; ++i)
        {
            items[i] = items[i + extraItems];
        }
    }

    public final List<BachelorGermlineVariant> getBachelorVariants() { return mGermlineVariants; }

}

