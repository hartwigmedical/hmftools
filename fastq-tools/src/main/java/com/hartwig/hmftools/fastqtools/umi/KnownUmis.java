package com.hartwig.hmftools.fastqtools.umi;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BASE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_DELIM;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_START;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_BASES;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_ID;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_QUALS;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;

public class KnownUmis
{
    private final String mUmiDelim;
    private final int mKnownUmiBaseDiff;
    private final boolean mKnownUmiUseNumeric;

    private final List<String> mKnownUmis;
    private final List<String> mKnownUmisReversed;
    private final AtomicLong mKnownUmiMatchCount;

    private final String UNMATCHED_UMI = "X";

    public KnownUmis(final UmiConfig config)
    {
        this(config.KnownUmiFile, config.UmiDelim, config.KnownUmiBaseDiff, config.KnownUmiUseNumeric);
    }

    public KnownUmis(final String knownUmiFile, final String umiDelim, final int knownUmiBaseDiff, final boolean knownUmiUseNumeric)
    {
        mUmiDelim = umiDelim;
        mKnownUmiBaseDiff = knownUmiBaseDiff;
        mKnownUmiUseNumeric = knownUmiUseNumeric;

        mKnownUmis = Lists.newArrayList();
        mKnownUmisReversed = Lists.newArrayList();
        mKnownUmiMatchCount = new AtomicLong();

        if(knownUmiFile != null)
        {
            loadKnownUmis(knownUmiFile);
        }
    }

    public boolean enabled() { return !mKnownUmis.isEmpty(); }

    public void adjustWithKnownUmi(
            final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        // find the longest matching UMI

        String umiBases1 = findKnownUmiMatch(r1ReadBuffer[READ_ITEM_BASES]);
        String umiBases2 = findKnownUmiMatch(r2ReadBuffer[READ_ITEM_BASES]);

        // make UMIs of the form 4+5, where integers are the length of the UMIs
        int umiLength1 = umiBases1.length();
        int umiLength2 = umiBases2.length();

        if(umiLength1 > 0 && umiLength2 > 0)
        {
            mKnownUmiMatchCount.incrementAndGet();
        }
        else
        {
            if(umiBases1.isEmpty())
            {
                umiBases1 = UNMATCHED_UMI;
                umiLength1 = 0;
            }

            if(umiBases2.isEmpty())
            {
                umiBases2 = UNMATCHED_UMI;
                umiLength2 = 0;
            }
        }

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId;

        if(mKnownUmiUseNumeric)
        {
            duplexUmiId = umiLength1 + mUmiDelim + umiLength2;
        }
        else
        {
            duplexUmiId = umiBases1 + mUmiDelim + umiBases2;
        }

        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(umiLength2);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(umiLength2);
    }

    private String findKnownUmiMatch(final String readBases)
    {
        for(int i = 0; i < mKnownUmis.size(); ++i)
        {
            String umi = mKnownUmis.get(i);
            String readUmi = readBases.substring(0, umi.length());

            if(readUmi.equals(umi))
                return umi;

            String umiReversed = mKnownUmisReversed.get(i);

            if(readUmi.equals(umiReversed))
                return umiReversed;
        }

        int knownUmiBaseDiff = mKnownUmiBaseDiff;

        if(readBases.charAt(0) == DNA_N_BASE)
            knownUmiBaseDiff = max(knownUmiBaseDiff, 1);

        if(knownUmiBaseDiff == 0)
            return "";

        // check for a 1-base mismatch
        for(int i = 0; i < mKnownUmis.size(); ++i)
        {
            String umi = mKnownUmis.get(i);
            String readUmi = readBases.substring(0, umi.length());

            if(!exceedsUmiDiff(readUmi, umi, knownUmiBaseDiff))
                return umi;

            String umiReversed = mKnownUmisReversed.get(i);

            if(!exceedsUmiDiff(readUmi, umiReversed, knownUmiBaseDiff))
                return umiReversed;
        }

        return "";
    }

    protected static boolean exceedsUmiDiff(final String first, final String second, int permittedDiff)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > permittedDiff)
                    return true;
            }
        }

        return false;
    }

    private void loadKnownUmis(final String knownUmiFile)
    {
        List<String> knownUmis = loadDelimitedIdFile(knownUmiFile, "KnownUmi", CSV_DELIM);
        knownUmis.forEach(x -> mKnownUmis.add(x));

        Collections.sort(mKnownUmis, Comparator.comparingInt(x -> -x.length()));

        mKnownUmis.forEach(x -> mKnownUmisReversed.add(Nucleotides.reverseComplementBases(x)));

        if(!knownUmis.isEmpty())
        {
            FQ_LOGGER.info("loaded {} known UMIs from {}", mKnownUmis.size(), knownUmiFile);
        }
    }

    public void logResults(final long readCount)
    {
        if(mKnownUmis.isEmpty() || readCount == 0)
            return;

        long matchCount = mKnownUmiMatchCount.get();
        double matchedPerc = matchCount / (double)readCount;
        FQ_LOGGER.info("known UMI matched({} {}%)", matchCount, format("%.2f", matchedPerc * 100));
    }
}
