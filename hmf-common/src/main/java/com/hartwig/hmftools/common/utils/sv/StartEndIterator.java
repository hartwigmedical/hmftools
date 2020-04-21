package com.hartwig.hmftools.common.utils.sv;

import java.util.List;

import com.google.common.collect.Lists;

public class StartEndIterator
{
    // iterators for start and end data
    public static final int SE_START = 0;
    public static final int SE_END = 1;
    public static final int SE_PAIR = 2;

    public static final String START_STR = "start";
    public static final String END_STR = "end";

    public static boolean isStart(int iter) { return iter == SE_START; }

    public static int seIndex(boolean isStart) { return isStart ? SE_START : SE_END; }

    public static int switchIndex(int iter) { return iter == SE_START ? SE_END : SE_START; }

}
