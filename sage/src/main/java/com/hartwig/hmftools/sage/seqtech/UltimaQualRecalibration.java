package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getColumnIndex;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class UltimaQualRecalibration
{
    // T0_RECAL	AC_T	56
    // TP_RECAL	1_A_FALSE	49

    public static final String CFG_FILENAME = "ultima_bqr_file";

    private final static byte RECALIBRATED_QUAL_MAX = 55;
    private final static int RECALIBRATED_QUAL_MAX_HP_LENGTH = 5;

    private final Map<String,Double> mTpQualMap;
    private final Map<String,Double> mT0QualMap;
    private byte mMaxRawQual;

    public UltimaQualRecalibration()
    {
        mTpQualMap = Maps.newHashMap();
        mT0QualMap = Maps.newHashMap();
        mMaxRawQual = 0;
    }

    public void setMaxRawQual(byte qual) { mMaxRawQual = qual; }
    public byte maxRawQual() { return mMaxRawQual; }

    public byte calcTpRecalibratedQual(final byte readQual, final int homoploymerLength, final char base, final boolean tpIsZero)
    {
        if(readQual < ULTIMA_MAX_QUAL || homoploymerLength > RECALIBRATED_QUAL_MAX_HP_LENGTH)
            return readQual;

        double bqrQual = getTpRecalibratedQual(homoploymerLength, base, tpIsZero);
        return bqrQual != INVALID_BASE_QUAL ? (byte)round(bqrQual) : readQual;
    }

    public double getTpRecalibratedQual(final int homoploymerLength, final char base, final boolean tpIsZero)
    {
        String key = TpData.formKey(homoploymerLength, base, tpIsZero);
        Double recalibratedQual = mTpQualMap.get(key);
        return recalibratedQual != null ? recalibratedQual : INVALID_BASE_QUAL;
    }

    public double getT0RecalibratedQual(final String triNucContext, final char base)
    {
        String key = T0Data.formKey(triNucContext, base);
        Double recalibratedQual = mT0QualMap.get(key);
        return recalibratedQual != null ? recalibratedQual : RECALIBRATED_QUAL_MAX;
    }

    private enum QualType { T0, TP; }

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
                double recalibratedQual = Double.parseDouble(values[rqIndex]);

                if(type == QualType.TP)
                {
                    mTpQualMap.put(key, recalibratedQual);
                }
                else
                {
                    mT0QualMap.put(key, recalibratedQual);
                }
            }
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to read Ultima base-qual recalibration file({}): {}", filename, e.toString());
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
            return format("%d_%c_%s", homoploymerLength, base, tpIsZero);
        }

        public String key() { return mKey; }

        public String toString() { return format("%d-%c %s %.0f", HomoploymerLength, Base, TpIsZero ? "tp-zero" : "tp-positive", RecalibratedQual); }
    }
}
