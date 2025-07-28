package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NeoCommon
{
    public static final String APP_NAME = "Neo";

    public static final int DOWNSTREAM_PRE_GENE_DISTANCE = 100000; // in concordance with Linx

    public static final Logger NE_LOGGER = LogManager.getLogger(NeoCommon.class);

    public static final List<String> IMMUNE_TRANSCRIPT_PREFIXES = Lists.newArrayList(
            "TR_J", "TR_V", "TR_D", "IG_J", "IG_D", "IG_J");

    public static String transcriptsToStr(final Collection<String> transcripts)
    {
        final StringJoiner transStr = new StringJoiner(ITEM_DELIM);
        transcripts.forEach(x -> transStr.add(x));
        return transStr.toString();
    }
}
