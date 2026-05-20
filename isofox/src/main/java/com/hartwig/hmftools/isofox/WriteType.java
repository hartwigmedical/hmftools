package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public enum WriteType
{
    EXON, // exon region data
    SPLICE_JUNC, // canonical splice junction counts
    FRAG_LENGTH, // intronic fragment lengths
    FRAG_LENGTH_BY_GENE, // intronic fragment lengths by gene
    READ, // all read attribution
    SPLICE_SITE, // splice site support
    TRANS_COMBO, // transcript group data for EM algo
    GC_RATIO; // GC ratio counts from all genic reads

    public static List<WriteType> parseConfig(final String configStr)
    {
        if(configStr == null || configStr.isEmpty())
            return Collections.emptyList();

        String[] configItems = configStr.split(ITEM_DELIM, -1);
        List<WriteType> writeTypes = Lists.newArrayList();

        for(String configItem : configItems)
        {
            writeTypes.add(WriteType.valueOf(configItem));
        }

        return writeTypes;
    }
}
