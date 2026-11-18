package com.hartwig.hmftools.taur.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.taur.FastaCommon.READ_ID_BREAK;
import static com.hartwig.hmftools.taur.FastaCommon.READ_ID_DELIM;
import static com.hartwig.hmftools.taur.FastaCommon.READ_ID_START;
import static com.hartwig.hmftools.taur.FastaCommon.READ_ITEM_BASES;
import static com.hartwig.hmftools.taur.FastaCommon.READ_ITEM_ID;
import static com.hartwig.hmftools.taur.FastaCommon.READ_ITEM_QUALS;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.taur.TaurConfig;

public final class UmiExtractor
{
    private final int mUmiLength;
    private final String mUmiDelim;
    private final int mAdapterLength;
    private final int mAdapterUmiLength;
    private final String mAdapterSequence;
    public final String mAdapterSequenceReversed;

    private final KnownUmis mKnownUmis;

    public UmiExtractor(
            final int umiLength, final String umiDelim, final int adapterLength, final String adapterSequence,
            final String knownUmiFile, final int knownUmiBaseDiff, final boolean knownUmiUseNumeric)
    {
        mUmiLength = umiLength;
        mUmiDelim = umiDelim;
        mAdapterLength = adapterLength;
        mAdapterSequence = adapterSequence;

        if(mAdapterSequence != null)
        {
            mAdapterUmiLength = mAdapterSequence.length() + mUmiLength;
            mAdapterSequenceReversed = Nucleotides.reverseComplementBases(mAdapterSequence);
        }
        else
        {
            mAdapterUmiLength = mAdapterLength > 0 ? mUmiLength + mAdapterLength : mUmiLength;
            mAdapterSequenceReversed = null;
        }

        mKnownUmis = new KnownUmis(knownUmiFile, umiLength, umiDelim, knownUmiBaseDiff, knownUmiUseNumeric);
    }

    public UmiExtractor(final TaurConfig config)
    {
        this(config.UmiLength, config.UmiDelim, config.AdapterLength, config.AdapterSequence,
                config.KnownUmiFile, config.KnownUmiBaseDiff, config.KnownUmiUseNumeric);
    }

    public boolean processReadBases(final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        /* expected format:
        @A00121:853:H5JNNMMABC:1:1101:1443:1047 1:N:0:GGCACAACCT+CAGGAGTCTA
        GNGAGATGGAGAATTTTCTGGAGATGTCTGAGGAATTTTTTCCTCAGTCTTAAGAGTA etc
        +
        F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF etc
        */

        // data validation
        if(r1ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START || r2ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START)
            return false;

        int minReadLength = max(mAdapterLength, mUmiLength);

        if(r1ReadBuffer[READ_ITEM_BASES].length() <= minReadLength || r1ReadBuffer[READ_ITEM_QUALS].length() <= minReadLength)
            return false;

        int delimIndex = r1ReadBuffer[READ_ITEM_ID].indexOf(READ_ID_BREAK);

        if(delimIndex < 1 || r2ReadBuffer[READ_ITEM_ID].charAt(delimIndex) != READ_ID_BREAK)
            return false;

        String readId1 = r1ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);
        String readId2 = r2ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);

        if(!readId1.equals(readId2))
            return false;

        if(mAdapterLength > 0)
        {
            adjustWithAdapter(readId1, readId2, delimIndex, r1ReadBuffer,  r2ReadBuffer);
        }
        else if(mKnownUmis.enabled())
        {
            mKnownUmis.adjustWithKnownUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
        }
        else
        {
            adjustWithFixedUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
        }

        return true;
    }

    public void adjustWithFixedUmi(final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        String umiBases1 = r1ReadBuffer[READ_ITEM_BASES].substring(0, mUmiLength);
        String umiBases2 = r2ReadBuffer[READ_ITEM_BASES].substring(0, mUmiLength);

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId = umiBases1 + mUmiDelim + umiBases2;
        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(mUmiLength);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(mUmiLength);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(mUmiLength);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(mUmiLength);
    }

    private void adjustWithAdapter(
            final String readId1, final String readId2, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        String adapterUmiBases1 = r1ReadBuffer[READ_ITEM_BASES].substring(0, mAdapterUmiLength);

        String umiBases1 = adapterUmiBases1.substring(0, mUmiLength);

        // append UMIs to read Id and remove from bases and quals
        String newReadId = readId1 + READ_ID_DELIM + umiBases1;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(mAdapterUmiLength);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(mAdapterUmiLength);

        // the R2 read may have the reversed adapter+UMI sequence at the end
        String adapterUmiBases2Start = r2ReadBuffer[READ_ITEM_BASES].substring(0, mAdapterUmiLength);
        int read2Length = r2ReadBuffer[READ_ITEM_BASES].length();
        String adapterUmiBases2End = r2ReadBuffer[READ_ITEM_BASES].substring(read2Length - mAdapterUmiLength);

        int adapterSeqIndex = adapterUmiBases2Start.indexOf(mAdapterSequence);

        if(adapterSeqIndex >= 0)
        {
            // trim from start
            int trimIndex = adapterSeqIndex + mAdapterSequence.length();
            r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(trimIndex);
            r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(trimIndex);
        }
        else
        {
            adapterSeqIndex = adapterUmiBases2End.indexOf(mAdapterSequenceReversed);

            if(adapterSeqIndex >= 0)
            {
                int trimIndex = read2Length - mAdapterUmiLength + adapterSeqIndex;
                r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(0, trimIndex);
                r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(0, trimIndex);
            }
        }
    }

    public void logResults(final long readCount)
    {
        mKnownUmis.logResults(readCount);
    }
}
