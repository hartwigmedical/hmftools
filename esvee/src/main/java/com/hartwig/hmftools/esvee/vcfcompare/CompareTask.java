package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

enum CompareTask
{
    MATCH_BREAKENDS,
    LINE_COMPARE;

    public static List<CompareTask> getAllTasks()
    {
        return List.of(CompareTask.class.getEnumConstants());
    }

    private static final String ALL = "ALL";

    public static List<CompareTask> fromConfig(String configStr)
    {
        if(configStr == null || configStr.equals(ALL))
            return getAllTasks();

        Set<CompareTask> compareTasks = new HashSet<>();

        String[] configStrValues = configStr.split(ITEM_DELIM, -1);

        for(String value : configStrValues)
            compareTasks.add(CompareTask.valueOf(value));

        return new ArrayList<>(compareTasks);
    }
}
