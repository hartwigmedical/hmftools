package com.hartwig.hmftools.orange.cohort.application;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Record;
import org.jooq.Result;

public final class SampleQuery
{
    private static final String DOID_SEPARATOR = ",";

    @NotNull
    public static List<Sample> selectFromDatarequest(@NotNull DatabaseAccess database)
    {
        return select(database, "datarequest");
    }

    @NotNull
    public static List<Sample> selectFromClinical(@NotNull DatabaseAccess database)
    {
        return select(database, "clinical");
    }

    @NotNull
    private static List<Sample> select(@NotNull DatabaseAccess database, @NotNull String table)
    {
        List<Sample> samples = Lists.newArrayList();

        Result<Record> result = database.context().resultQuery("select sampleId, doids from " + table).fetch();
        for(Record record : result)
        {
            String sampleId = (String) record.getValue(0);
            String doidString = (String) record.getValue(1);
            if(doidString == null)
            {
                LOGGER.debug(" Skipping sample {} because doids are unknown", sampleId);
            }
            else
            {
                samples.add(ImmutableSample.builder().sampleId(sampleId).doids(toDoids(doidString)).build());
            }
        }

        return samples;
    }

    @NotNull
    private static Set<String> toDoids(@Nullable String doidString)
    {
        if(doidString.isEmpty())
        {
            return Sets.newHashSet();
        }

        return Sets.newHashSet(doidString.split(DOID_SEPARATOR));
    }
}
