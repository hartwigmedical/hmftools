package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getColumnIndex;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMRecord;

public class UltimaQualRecalibration
{
    // T0_RECAL	AC_T	56
    // TP_RECAL	1_A_FALSE	49

    public static final String CFG_FILENAME = "ultima_bqr_file";

    private final static String T0_OUT_OF_CYCLE = "OUT_OF_CYCLE";
    private final static int RECALIBRATED_QUAL_MAX_HP_LENGTH = 5;
    protected final static byte MAX_TP_T0_RECALIBRATION = 50;

    private final Map<String,Byte> mTpQualMap;
    private final Map<String,Byte> mT0QualMap;
    private byte mMaxRawQual;
    private byte mOutOfCycleT0Qual;

    public UltimaQualRecalibration()
    {
        mTpQualMap = Maps.newHashMap();
        mT0QualMap = Maps.newHashMap();
        mMaxRawQual = ULTIMA_MAX_QUAL;
        mOutOfCycleT0Qual = MAX_TP_T0_RECALIBRATION;
    }

    public void setMaxRawQual(byte qual) { mMaxRawQual = qual; }
    public byte maxRawQual() { return mMaxRawQual; }
    public byte outOfCycleT0Qual() { return mOutOfCycleT0Qual; }

    public byte calcTpRecalibratedQual(
            final byte readQual, final int homoploymerLength, final char base, final boolean tpIsZero, final boolean isReverse)
    {
        if(readQual < mMaxRawQual || homoploymerLength > RECALIBRATED_QUAL_MAX_HP_LENGTH)
            return readQual;

        byte bqrQual = getTpRecalibratedQual(homoploymerLength, base, tpIsZero, isReverse);
        return bqrQual != INVALID_BASE_QUAL ? bqrQual : readQual;
    }

    public byte getTpRecalibratedQual(final int homoploymerLength, final char base, final boolean tpIsZero, final boolean isReverse)
    {
        char baseToLookup = isReverse ? swapDnaBase(base) : base;
        String key = TpData.formKey(homoploymerLength, baseToLookup, tpIsZero);
        Byte recalibratedQual = mTpQualMap.get(key);
        return recalibratedQual != null ? recalibratedQual : INVALID_BASE_QUAL;
    }

    public byte calcT0RecalibratedQual(final byte readQual, final SAMRecord record, int varReadIndex)
    {
        if(readQual < mMaxRawQual)
            return readQual;

        byte bqrQual = getT0RecalibratedQual(record, varReadIndex);
        return bqrQual != INVALID_BASE_QUAL ? bqrQual : readQual;
    }

    public byte getT0RecalibratedQual(final SAMRecord record, int varReadIndex)
    {
        if(varReadIndex <= 0 || varReadIndex >= record.getReadBases().length - 1)
            return ULTIMA_INVALID_QUAL;

        char variantBase = (char)record.getReadBases()[varReadIndex];
        boolean isReverse = record.getReadNegativeStrandFlag();
        char baseToLookup = isReverse ? swapDnaBase(variantBase) : variantBase;
        String tnc = (char)record.getReadBases()[varReadIndex - 1] + String.valueOf((char)record.getReadBases()[varReadIndex + 1]);
        String tncToLookup = isReverse ? reverseComplementBases(tnc) : tnc;
        return getT0RecalibratedQual(tncToLookup, baseToLookup);
    }

    public byte getT0RecalibratedQual(final String triNucContext, final char base)
    {
        String key = T0Data.formKey(triNucContext, base);
        Byte recalibratedQual = mT0QualMap.get(key);
        return recalibratedQual != null ? recalibratedQual : mOutOfCycleT0Qual;
    }

    private enum QualType { T0_RECAL, TP_RECAL; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CFG_FILENAME, false, "Ultima base qual recalibration file");
    }

    private static final String COL_TYPE = "Type";
    private static final String COL_KEY = "Key";
    private static final String COL_RECAL_QUAL = "RecalibratedQual";

    public void loadRecalibrationFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());
            loadRecalibrationData(lines);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to read Ultima base-qual recalibration file({}): {}", filename, e.toString());
        }
    }

    protected void loadRecalibrationData(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        int typeIndex = getColumnIndex(fieldsIndexMap, COL_TYPE);
        int keyIndex = getColumnIndex(fieldsIndexMap, COL_KEY);
        int rqIndex = getColumnIndex(fieldsIndexMap, COL_RECAL_QUAL);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            QualType type = QualType.valueOf(values[typeIndex]);
            String key = values[keyIndex];
            byte recalibratedQual = (byte) Math.min(Byte.parseByte(values[rqIndex]), MAX_TP_T0_RECALIBRATION);

            if(type == QualType.TP_RECAL)
            {
                mTpQualMap.put(key, recalibratedQual);
            }
            else
            {
                if(key.equals(T0_OUT_OF_CYCLE))
                    mOutOfCycleT0Qual = recalibratedQual;
                else
                    mT0QualMap.put(key, recalibratedQual);
            }
        }
    }

    private class T0Data
    {
        public final String TriNucContext;
        public final char Base;
        public final double RecalibratedQual;

        private final String mKey;

        public T0Data(final String triNucContext, final char base, final double recalibratedQual)
        {
            TriNucContext = triNucContext;
            Base = base;
            RecalibratedQual = recalibratedQual;
            mKey = formKey(TriNucContext, Base);
        }

        public String key() { return mKey; }

        public static String formKey(final String triNucContext, final char base)
        {
            return format("%s_%c", triNucContext, base);
        }

        public String toString() { return format("%s=%.0f", mKey, RecalibratedQual); }
    }

    private class TpData
    {
        public final int HomoploymerLength;
        public final char Base;
        public final boolean TpIsZero;
        public final double RecalibratedQual;

        private final String mKey;

        public TpData(final int homoploymerLength, final char base, final boolean tpIsZero, final double recalibratedQual)
        {
            HomoploymerLength = homoploymerLength;
            Base = base;
            TpIsZero = tpIsZero;
            RecalibratedQual = recalibratedQual;
            mKey = formKey(HomoploymerLength, Base, TpIsZero);
        }

        public static String formKey(final int homoploymerLength, final char base, final boolean tpIsZero)
        {
            return format("%d_%c_%s", homoploymerLength, base, String.valueOf(tpIsZero).toUpperCase());
        }

        public String key() { return mKey; }

        public String toString() { return format("%d-%c %s %.0f", HomoploymerLength, Base, TpIsZero ? "tp-zero" : "tp-positive", RecalibratedQual); }
    }
}
