package com.hartwig.hmftools.fastqtools.umi;

import static java.lang.String.format;

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

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;

public class KnownUmis
{
    private final UmiConfig mConfig;
    private final List<String> mKnownUmis;
    private final List<String> mKnownUmisReversed;
    private long mKnownUmiMatchCount;

    public KnownUmis(final UmiConfig config)
    {
        mConfig = config;

        mKnownUmis = Lists.newArrayList();
        mKnownUmisReversed = Lists.newArrayList();
        mKnownUmiMatchCount = 0;

        if(mConfig.KnownUmiFile != null)
        {
            loadKnownUmis(mConfig.KnownUmiFile);
        }
    }

    public boolean enabled() { return !mKnownUmis.isEmpty(); }

    public void adjustWithKnownUmi(
            final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        // find the longest matching UMI

        String umiBases1 = findKnownUmiMatch(r1ReadBuffer[READ_ITEM_BASES], false);
        String umiBases2 = findKnownUmiMatch(r2ReadBuffer[READ_ITEM_BASES], true);

        // make UMIs of the form 4+5, where integers are the length of the UMIs
        int umiLength1 = umiBases1.length();
        int umiLength2 = umiBases2.length();

        if(umiLength1 > 0 && umiLength2 > 0)
            ++mKnownUmiMatchCount;

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId = umiLength1 + mConfig.UmiDelim + umiLength2;
        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(umiLength2);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(umiLength2);
    }

    private String findKnownUmiMatch(final String readBases, boolean requiresReverse)
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

        if(mConfig.KnownUmiBaseDiff == 0)
            return "";

        // check for a 1-base mismatch
        for(int i = 0; i < mKnownUmis.size(); ++i)
        {
            String umi = mKnownUmis.get(i);
            String readUmi = readBases.substring(0, umi.length());

            if(!exceedsUmiDiff(readUmi, umi, mConfig.KnownUmiBaseDiff))
                return umi;

            String umiReversed = mKnownUmisReversed.get(i);

            if(!exceedsUmiDiff(readUmi, umiReversed, mConfig.KnownUmiBaseDiff))
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

        double matchedPerc = mKnownUmiMatchCount / (double)readCount;
        FQ_LOGGER.info("known UMI matched({} {}%)", mKnownUmiMatchCount, format("%.2f", matchedPerc * 100));
    }
}
