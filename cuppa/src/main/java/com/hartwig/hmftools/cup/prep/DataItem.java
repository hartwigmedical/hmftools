package com.hartwig.hmftools.cup.prep;

import java.util.Comparator;
import java.util.Objects;

public class DataItem
{
    public final Index Index;
    public final String Value;

    public static final String FLD_SOURCE = "Source";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_KEY = "Key";
    public static final String FLD_VALUE = "Value";

    public static class Index
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
    }

    public DataItem(final DataSource source, final ItemType type, final String key, final String value)
    {
        Index = new Index(source, type, key);
        Value = value;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;

        if(o == null || getClass() != o.getClass())
            return false;

        final DataItem dataItem = (DataItem) o;
        return Index.Source == dataItem.Index.Source
                && Index.Type == dataItem.Index.Type
                && Index.Key.equals(dataItem.Index.Key)
                && Value.equals(dataItem.Value);
    }

    @Override
    public String toString()
    {
        return Index + ", Value=" + Value;
    }

    public static class IndexComparator implements Comparator<Index>
    {
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
        public int compare(Index index1, Index index2)
        {
            int typeDiff = index1.Type.compareTo(index2.Type);
            if(typeDiff != 0)
                return typeDiff;

            if(index1.Type == ItemType.GEN_POS)
            {
                return compareGenomicLoci(index1.Key, index2.Key, "_");
            }

            if(index1.Type == ItemType.ALT_SJ)
            {
                return compareGenomicLoci(index1.Key, index2.Key, ";");
            }

            if(index1.Type == ItemType.SIGNATURE)
            {
                int sigNum1 = Integer.parseInt(index1.Key.split("_")[1]);
                int sigNum2 = Integer.parseInt(index2.Key.split("_")[1]);
                return sigNum1 - sigNum2;
            }

            return index1.Key.compareTo(index2.Key);
        }
    }
}
