package com.hartwig.hmftools.taur;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class FastaCommon
{
    public static final Logger TR_LOGGER = LogManager.getLogger(FastaCommon.class);

    public static final int READ_ITEM_ID = 0;
    public static final int READ_ITEM_BASES = 1;
    public static final int READ_ITEM_SPARE = 2;
    public static final int READ_ITEM_QUALS = 3;
    public static final int READ_LINE_COUNT = 4;

    public static final char READ_ID_START = '@';
    public static final char READ_ID_BREAK = ' ';
    public static final char READ_ID_DELIM = ':';

    public static final String FASTQ_SUFFIX_STANDARD = "fastq";
    public static final String FASTQ_ZIP_EXTENSION = ".fastq.gz";
    public static final String FASTQ_SUFFIX_SHORT = "fq";
}
