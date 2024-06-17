package com.hartwig.hmftools.cup.common;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CupConstants
{
    public static final String APP_NAME = "Cuppa";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CupConstants.class);

    public static final int GEN_POS_BUCKET_SIZE = 500000;

    public static final List<String> AID_APOBEC_TRINUCLEOTIDE_CONTEXTS = Lists.newArrayList(
            "C>T_TCA", "C>T_TCC", "C>T_TCG", "C>T_TCT", "C>G_TCA", "C>G_TCC", "C>G_TCG", "C>G_TCT");
}
