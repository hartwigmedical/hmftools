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

    public class Index
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
            {
                return true;
            }
            if(o == null || getClass() != o.getClass())
            {
                return false;
            }
            final Index index = (Index) o;
            return Source == index.Source && Type == index.Type && Key.equals(index.Key);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(Source, Type, Key);
        }

        public String toString()
        {
            return String.format(
                    "%s=%s, %s=%s, %s=%s",
                    FLD_SOURCE, Index.Source,
                    FLD_CATEGORY, Index.Type,
                    FLD_KEY, Index.Key
            );
        }
    }

    public DataItem(final DataSource source, final ItemType type, final String key, final String value)
    {
        Index = new Index(source, type, key);
        Value = value;
    }

    public static class IndexComparator implements Comparator<Index>
    {
        private String padChrom(String chrom)
        {
//            chrom = chrom.replace("chr", "");
            char padChar = (Character.isDigit(chrom.charAt(0))) ? '0' : 'a';
            return String.format("%2s", chrom).replace(' ', padChar);
        }

        @Override
        public int compare(Index index1, Index index2)
        {
            int typeDiff = index1.Type.compareTo(index2.Type);
            if(typeDiff != 0)
                return typeDiff;

            if(index1.Type == ItemType.GEN_POS)
            {
                String[] key1 = index1.Key.split("_");
                String[] key2 = index2.Key.split("_");

                String chrom1 = padChrom(key1[0]);
                String chrom2 = padChrom(key2[0]);

                int chromDiff = chrom1.compareTo(chrom2);
                if(chromDiff != 0)
                    return chromDiff;

                Integer pos1 = Integer.parseInt(key1[1]);
                Integer pos2 = Integer.parseInt(key2[1]);

                return Integer.signum(pos1 - pos2);
            }

            return index1.Key.compareTo(index2.Key);
        }
    }
}
