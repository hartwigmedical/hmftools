package com.hartwig.hmftools.compar.resource;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BED_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ZIP_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.compar.resource.ResourceDefinition.BEDPE_COLUMNS;
import static com.hartwig.hmftools.compar.resource.ResourceDefinition.BED_COLUMNS;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.compar.common.field.FieldValue;

public class ResourceFile
{
    private final String mDelim;
    private final List<String> mFileColumns;
    private final int mInfoIndex;
    private final Set<String> mInfoKeys;
    private final Map<String,List<String>> mLinesByKey;
    private final List<ResourceField> mFields;
    private final int[] mFieldIndices;

    private int mLineNumber;

    private static final String META_PREFIX = "##";
    private static final String HEADER_PREFIX = "#";
    private static final String BEDPE_EXTENSION = ".bedpe";

    private static final String INFO_DELIM = ";";
    private static final String INFO_MISSING = ".";
    private static final char INFO_VALUE_DELIM = '=';

    public ResourceFile(final ResourceDefinition definition, final Path path) throws IOException
    {
        String fileName = path.getFileName().toString();
        mDelim = inferFileDelimiter(fileName);
        mInfoKeys = Sets.newTreeSet();
        mLinesByKey = Maps.newLinkedHashMap();
        mLineNumber = 0;

        try(BufferedReader reader = createBufferedReader(path.toString()))
        {
            if(fileName.endsWith(BEDPE_EXTENSION) || fileName.endsWith(BEDPE_EXTENSION + ZIP_EXTENSION))
                mFileColumns = BEDPE_COLUMNS;
            else if(fileName.endsWith(BED_EXTENSION) || fileName.endsWith(BED_EXTENSION + ZIP_EXTENSION))
                mFileColumns = BED_COLUMNS;
            else
                mFileColumns = readHeader(reader);

            mInfoIndex = definition.InfoColumn != null ? mFileColumns.indexOf(definition.InfoColumn) : -1;
            checkKeyColumns(definition);
            readRows(definition, reader);
        }

        List<String> columns = columns();
        mFields = definition.Fields.stream().filter(x -> x.DefaultValue != null || columns.contains(x.name())).collect(Collectors.toList());
        mFieldIndices = mFields.stream().mapToInt(x -> mFileColumns.indexOf(x.name())).toArray();
    }

    public Map<String,List<String>> linesByKey() { return mLinesByKey; }

    public List<String> fileColumns() { return mFileColumns; }

    public List<String> columns()
    {
        if(mInfoIndex < 0)
            return mFileColumns;

        List<String> columns = Lists.newArrayList(mFileColumns);
        columns.remove(mInfoIndex);
        columns.addAll(mInfoKeys);
        return columns;
    }

    public Map<String,FieldValue> fieldValues(final String line)
    {
        String[] values = line.split(mDelim, mFileColumns.size());
        Map<String,String> info = mInfoIndex >= 0 ? parseInfo(values[mInfoIndex]) : Collections.emptyMap();

        Map<String,FieldValue> fieldValues = Maps.newLinkedHashMap();

        for(int i = 0; i < mFields.size(); ++i)
        {
            ResourceField field = mFields.get(i);
            String value;

            if(mFieldIndices[i] >= 0)
                value = values[mFieldIndices[i]];
            else if(field.DefaultValue != null)
                value = field.DefaultValue;
            else
                value = info.getOrDefault(field.name(), "");

            fieldValues.put(field.name(), field.createValue(value));
        }

        return fieldValues;
    }

    private List<String> readHeader(final BufferedReader reader) throws IOException
    {
        String header = readLine(reader);

        while(header != null && header.startsWith(META_PREFIX))
        {
            header = readLine(reader);
        }

        if(header == null)
            throw new IOException("missing header");

        if(header.startsWith(HEADER_PREFIX))
            header = header.substring(HEADER_PREFIX.length());

        return Arrays.asList(header.split(mDelim, -1));
    }

    private void checkKeyColumns(final ResourceDefinition definition) throws IOException
    {
        List<String> missingColumns = definition.KeyColumns.stream().filter(x -> !mFileColumns.contains(x)).collect(Collectors.toList());

        if(!missingColumns.isEmpty())
            throw new IOException(format("missing columns(%s)", String.join(ITEM_DELIM, missingColumns)));
    }

    private void readRows(final ResourceDefinition definition, final BufferedReader reader) throws IOException
    {
        int[] keyIndices = definition.KeyColumns.stream().mapToInt(mFileColumns::indexOf).toArray();
        int stripChrPrefixIndex = definition.StripChrPrefixColumn != null ? mFileColumns.indexOf(definition.StripChrPrefixColumn) : -1;

        String line;

        while((line = readLine(reader)) != null)
        {
            String[] values = line.split(mDelim, mFileColumns.size());

            if(values.length < mFileColumns.size())
                throw new IOException(format("line(%d) has %d of %d columns", mLineNumber, values.length, mFileColumns.size()));

            if(mInfoIndex >= 0)
                mInfoKeys.addAll(parseInfo(values[mInfoIndex]).keySet());

            mLinesByKey.computeIfAbsent(key(values, keyIndices, stripChrPrefixIndex), k -> Lists.newArrayList()).add(line);
        }
    }

    private String readLine(final BufferedReader reader) throws IOException
    {
        ++mLineNumber;
        return reader.readLine();
    }

    private static String key(final String[] values, final int[] keyIndices, final int stripChrPrefixIndex)
    {
        StringJoiner key = new StringJoiner(ITEM_DELIM);

        for(int keyIndex : keyIndices)
        {
            key.add(keyIndex == stripChrPrefixIndex ? stripChrPrefix(values[keyIndex]) : values[keyIndex]);
        }

        return key.toString();
    }

    private static Map<String,String> parseInfo(final String info)
    {
        Map<String,String> values = Maps.newHashMap();

        if(info.isEmpty() || info.equals(INFO_MISSING))
            return values;

        for(String item : info.split(INFO_DELIM))
        {
            int index = item.indexOf(INFO_VALUE_DELIM);

            if(index < 0)
                values.put(item, Boolean.TRUE.toString());
            else
                values.put(item.substring(0, index), item.substring(index + 1));
        }

        return values;
    }
}
