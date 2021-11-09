package com.hartwig.hmftools.orange.cohort.application;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Record;
import org.jooq.Result;

public final class SampleDataQuery {

    private static final Logger LOGGER = LogManager.getLogger(SampleDataQuery.class);
    private static final String DOID_SEPARATOR = ",";

    private SampleDataQuery() {
    }

    @NotNull
    public static List<Sample> run(@NotNull DatabaseAccess database) {
        List<Sample> samples = Lists.newArrayList();

        Result<Record> result = database.context().resultQuery("select sampleId, doids from datarequest").fetch();
        for (Record record : result) {
            String sampleId = (String) record.getValue(0);
            String doidString = (String) record.getValue(1);
            if (doidString == null) {
                LOGGER.info(" Skipping sample {} because doids are unknown", sampleId);
            } else {
                samples.add(ImmutableSample.builder().sampleId(sampleId).doids(toDoids(doidString)).build());
            }
        }

        return samples;
    }

    @NotNull
    private static Set<String> toDoids(@Nullable String doidString) {
        if (doidString.isEmpty()) {
            return Sets.newHashSet();
        }

        return Sets.newHashSet(doidString.split(DOID_SEPARATOR));
    }
}
