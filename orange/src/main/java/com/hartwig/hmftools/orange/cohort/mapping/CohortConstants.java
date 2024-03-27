package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public final class CohortConstants
{
    public static final String COHORT_PAN_CANCER = "Pan-cancer";
    public static final String COHORT_OTHER = "Other";
    public static final String COHORT_UNKNOWN = "Unknown";
    public static final String COHORT_ESOPHAGUS = "Esophagus";
    public static final String COHORT_STOMACH = "Stomach";

    public static final List<Set<String>> DOID_COMBINATIONS_TO_MAP_TO_OTHER = Lists.newArrayList();
    public static final List<Set<String>> DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS = Lists.newArrayList();
    public static final List<Set<String>> DOID_COMBINATIONS_TO_MAP_TO_STOMACH = Lists.newArrayList();
    public static final Map<Set<String>, String> DOID_COMBINATION_MAP = new HashMap<>();

    static
    {
        // The combination of liver cancer and bile duct is occasionally used but not easily mapped.
        DOID_COMBINATIONS_TO_MAP_TO_OTHER.add(Sets.newHashSet("686", "4947"));

        // Kidney cancer is a tough one to spread out without affecting either Kidney or Urothelial tract itself.
        DOID_COMBINATIONS_TO_MAP_TO_OTHER.add(Sets.newHashSet("263"));

        // Combination of urethra cancer and renal cell cancer cannot easily be mapped.
        DOID_COMBINATIONS_TO_MAP_TO_OTHER.add(Sets.newHashSet("734", "2671", "4450"));

        // Combination of esophagus or stomach cancer & gastroesophageal system cancer should resolve to esophagus or stomach.
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080374", "5041"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080374", "1107"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080374", "4914"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080374", "3748"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080375", "5041"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080375", "1107"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("0080375", "4914"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("4944", "5041"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("4944", "1107"));
        DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS.add(Sets.newHashSet("4944", "4914"));

        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080374", "10534"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080374", "5517"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080374", "5516"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080374", "3717"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080375", "10534"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080375", "5517"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("0080375", "3717"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("4944", "10534"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("4944", "5517"));
        DOID_COMBINATIONS_TO_MAP_TO_STOMACH.add(Sets.newHashSet("4944", "3717"));

        for (Set<String> combination : DOID_COMBINATIONS_TO_MAP_TO_OTHER)
        {
            DOID_COMBINATION_MAP.put(combination, COHORT_OTHER);
        }

        for (Set<String> combination : DOID_COMBINATIONS_TO_MAP_TO_ESOPHAGUS)
        {
            DOID_COMBINATION_MAP.put(combination, COHORT_ESOPHAGUS);
        }

        for (Set<String> combination : DOID_COMBINATIONS_TO_MAP_TO_STOMACH)
        {
            DOID_COMBINATION_MAP.put(combination, COHORT_STOMACH);
        }
    }
}
