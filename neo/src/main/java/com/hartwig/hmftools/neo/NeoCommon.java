package com.hartwig.hmftools.neo;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NeoCommon
{
    public static final int DOWNSTREAM_PRE_GENE_DISTANCE = 100000; // in concordance with Linx

    public static final Logger NE_LOGGER = LogManager.getLogger(NeoCommon.class);

    public static final List<String> IMMUNE_TRANSCRIPT_PREFIXES = Lists.newArrayList(
            "TR_J", "TR_V", "TR_D", "IG_J", "IG_D", "IG_J");
}
