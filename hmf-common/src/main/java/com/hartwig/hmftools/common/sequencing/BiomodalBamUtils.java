package com.hartwig.hmftools.common.sequencing;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.BASE_MODIFICATIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import java.util.List;
import java.util.SortedSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class BiomodalBamUtils
{
    public static String MODC_ANNOTATION = "MODC";
    public static String MM_PREFIX = "C+C.";
    public static String MM_SUFFIX = ";";

    public static String getMMValueFromModCReadIndices(
            final byte[] readBases, @Nullable final SortedSet<Integer> modCReadIndices, boolean isForward)
    {
        if(modCReadIndices == null || modCReadIndices.isEmpty())
            return MM_PREFIX + MM_SUFFIX;

        byte targetBase = isForward ? (byte) 'C' : (byte) swapDnaBase('C');
        StringBuilder skipsStr = new StringBuilder();
        int skip = 0;
        int inc = isForward ? 1 : -1;
        int readIndex = isForward ? 0 : readBases.length - 1;
        for(; readIndex >= 0 && readIndex < readBases.length; readIndex += inc)
        {
            byte base = readBases[readIndex];
            if(base != targetBase)
                continue;

            boolean isModC = modCReadIndices.contains(readIndex);
            if(!isModC)
            {
                skip++;
                continue;
            }

            skipsStr.append(',');
            skipsStr.append(skip);
            skip = 0;
        }

        return MM_PREFIX + skipsStr + MM_SUFFIX;
    }

    public static SortedSet<Integer> getModCReadIndicesFromMMSkipValues(
            final byte[] readBases, @Nullable final List<Integer> skipValues, boolean isForward)
    {
        if(skipValues == null || skipValues.isEmpty())
            return Sets.newTreeSet();

        byte targetBase = isForward ? (byte) 'C' : (byte) swapDnaBase('C');
        SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
        int inc = isForward ? 1 : -1;
        int readIndex = isForward ? 0 : readBases.length - 1;
        for(int skip : skipValues)
        {
           int skipCount = 0;
           for(; readIndex >= 0 && readIndex < readBases.length; readIndex += inc)
           {
               byte base = readBases[readIndex];
               if(base != targetBase)
                   continue;

               if(skipCount == skip)
               {
                   modCReadIndices.add(readIndex);
                   readIndex += inc;
                   break;
               }

               skipCount++;
           }
        }

        return modCReadIndices;
    }

    public static List<Integer> getMMSkipValuesFromMMValue(final String value)
    {
        if(value.equals(MM_PREFIX + MM_SUFFIX))
            return Lists.newArrayList();

        List<Integer> skipValues = Lists.newArrayList();
        String skipValuesStr = value.substring(MM_PREFIX.length() + 1, value.length() - MM_SUFFIX.length());
        String[] components = skipValuesStr.split(",");
        for(String skipStr : components)
            skipValues.add(Integer.parseInt(skipStr));

        return skipValues;
    }

    public static SortedSet<Integer> getModCReadIndicesFromMMValue(final byte[] readBases, final String mmValue, boolean isForward)
    {
        List<Integer> skipValues = getMMSkipValuesFromMMValue(mmValue);
        return getModCReadIndicesFromMMSkipValues(readBases, skipValues, isForward);
    }

    public static SortedSet<Integer> getModCReadIndices(final SAMRecord read)
    {
        String mmValue = read.getStringAttribute(BASE_MODIFICATIONS_ATTRIBUTE);
        if(mmValue == null)
            throw new IllegalArgumentException(format("readName(%s) must have an MM attribute: %s", read.getReadName(), read.getSAMString()));

        return getModCReadIndicesFromMMValue(read.getReadBases(), mmValue, !read.getReadNegativeStrandFlag());
    }
}
