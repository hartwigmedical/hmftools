package com.hartwig.hmftools.cup.prep;

import java.util.Locale;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class DataItem implements Comparable<DataItem>
{
    public final Index Index;
    public final String Value;

    public static final String FLD_SOURCE = "Source";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_KEY = "Key";
    public static final String FLD_VALUE = "Value";

    public DataItem(final DataSource source, final ItemType type, final String key, final int intValue)
    {
        Index = new Index(source, type, key);
        Value = String.valueOf(intValue);
    }

    public DataItem(final DataSource source, final ItemType type, final String key, final boolean boolValue)
    {
        Index = new Index(source, type, key);
        Value = boolValue ? "1" : "0";
    }

    public DataItem(final DataSource source, final ItemType type, final String key, final double doubleValue, String numberFormat)
    {
        Index = new Index(source, type, key);
        Value = String.format(Locale.ENGLISH, numberFormat, doubleValue);
    }

    @VisibleForTesting
    public DataItem(final DataSource source, final ItemType type, final String key, final String stringValue)
    {
        Index = new Index(source, type, key);
        Value = stringValue;
    }

    public static class Index implements Comparable<Index>
    {
        public final DataSource Source;
        public final ItemType Type;
        public final String Key;

        public Index(final DataSource source, final ItemType type, final String key)
        {
            Source = source;
            Type = type;
            Key = key;
        }

        @Override
        public String toString()
        {
            return String.format(
                    "%s=%s, %s=%s, %s=%s",
                    FLD_SOURCE, Source,
                    FLD_CATEGORY, Type,
                    FLD_KEY, Key
            );
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;

            if(o == null || getClass() != o.getClass())
                return false;

            final Index index = (Index) o;
            return Source == index.Source && Type == index.Type && Key.equals(index.Key);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(Source, Type, Key);
        }

        private String padChrom(String chrom)
        {
            char padChar = (Character.isDigit(chrom.charAt(0))) ? '0' : 'a';
            return String.format("%2s", chrom).replace(' ', padChar);
        }

        private int compareGenomicLoci(String locus1, String locus2, String delimiter)
        {
            String[] strings1 = locus1.split(delimiter);
            String[] strings2 = locus2.split(delimiter);

            String chrom1 = padChrom(strings1[0]);
            String chrom2 = padChrom(strings2[0]);
            int chromDiff = chrom1.compareTo(chrom2);
            if(chromDiff != 0)
                return chromDiff;

            Integer startPos1 = Integer.parseInt(strings1[1]);
            Integer startPos2 = Integer.parseInt(strings2[1]);
            int startPosDiff = startPos1 - startPos2;
            if(startPosDiff != 0)
                return startPosDiff;

            if(strings1.length==3)
            {
                Integer endPos1 = Integer.parseInt(strings1[2]);
                Integer endPos2 = Integer.parseInt(strings2[2]);
                int endPosDiff = endPos1 - endPos2;
                if(endPosDiff != 0)
                    return endPosDiff;
            }

            return 0;
        }

        @Override
        public int compareTo(@NotNull final DataItem.Index otherIndex)
        {
            int typeDiff = this.Type.compareTo(otherIndex.Type);
            if(typeDiff != 0)
                return typeDiff;

            if(otherIndex.Type == ItemType.GEN_POS)
            {
                return compareGenomicLoci(this.Key, otherIndex.Key, "_");
            }

            if(this.Type == ItemType.ALT_SJ)
            {
                return compareGenomicLoci(this.Key, otherIndex.Key, ";");
            }

            if(this.Type == ItemType.SIGNATURE)
            {
                int sigNum1 = Integer.parseInt(this.Key.split("_")[1]);
                int sigNum2 = Integer.parseInt(otherIndex.Key.split("_")[1]);
                return sigNum1 - sigNum2;
            }

            return this.Key.compareTo(otherIndex.Key);
        }
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;

        if(o == null || getClass() != o.getClass())
            return false;

        final DataItem otherDataItem = (DataItem) o;
        return Index.equals(otherDataItem.Index) && Value.equals(otherDataItem.Value);
    }

    @Override
    public int compareTo(@NotNull final DataItem otherDataItem)
    {
        return this.Index.compareTo(otherDataItem.Index);
    }

    @Override
    public String toString()
    {
        return Index + ", Value=" + Value;
    }
}
