package com.hartwig.hmftools.sage.quality;

public class QualityRecalibrationRecord implements Comparable<QualityRecalibrationRecord>
{
    public final BaseQualityKey Key;

    public final int Count;
    public final double RecalibratedQuality;

    public QualityRecalibrationRecord(final BaseQualityKey key, final int count, final double recalibratedQual)
    {
        Key = key;
        Count = count;
        RecalibratedQuality = recalibratedQual;
    }

    @Override
    public int compareTo(final QualityRecalibrationRecord other)
    {
        int countCompare = Integer.compare(other.Count, Count);
        if(countCompare != 0)
            return countCompare;

        int refCompare = Byte.compare(Key.Ref, other.Key.Ref);
        if(refCompare != 0)
            return refCompare;

        int altCompare = Byte.compare(Key.Alt, other.Key.Alt);
        if(altCompare != 0)
            return altCompare;

        int qualCompare = Byte.compare(Key.Quality, other.Key.Quality);
        if(qualCompare != 0)
            return qualCompare;

        byte[] context1 = Key.TrinucleotideContext;
        byte[] context2 = other.Key.TrinucleotideContext;

        if(context1 == null || context1.length < 3 || context2 == null || context2.length < 3)
        {
            return 0;
        }

        int triOne = Byte.compare(context1[0], context2[0]);
        if(triOne != 0)
            return triOne;

        int triTwo = Byte.compare(context1[1], context2[1]);
        if(triTwo != 0)
            return triTwo;

        int triThree = Byte.compare(context1[2], context2[2]);
        if(triThree != 0)
            return triThree;

        return 0;

    }

}
