package com.hartwig.hmftools.fastqtools.umi;

import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_DELIM;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_START;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_BASES;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_ID;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_QUALS;

import com.hartwig.hmftools.common.codon.Nucleotides;

public final class UmiExtractor
{
    private final int mUmiLength;
    private final String mUmiDelim;
    private final int mAdapterLength;
    private final int mAdapterUmiLength;
    private final String mAdapterSequence;
    public final String mAdapterSequenceReversed;

    public UmiExtractor(
            final int umiLength, final String umiDelim, final int adapterLength, final String adapterSequence)
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
    }

    public UmiExtractor(final UmiConfig config)
    {
        this(config.UmiLength, config.UmiDelim, config.AdapterLength, config.AdapterSequence);
    }

    public int adapterUmiLength() { return mAdapterUmiLength; }

    public void adjustWithFixedUmi(
            final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
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

    protected void adjustWithAdapter(
            final String readId1, final String readId2, final int delimIndex,
            final String[] r1ReadBuffer, final String[] r2ReadBuffer)
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
}
