package com.hartwig.hmftools.orange.cohort.application;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableObservation;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Record;
import org.jooq.Result;

public final class SvTmbQuery {

    private SvTmbQuery() {
    }

    @NotNull
    public static List<Observation> run(@NotNull DatabaseAccess database, @NotNull List<Sample> samplesToInclude) {
        List<Observation> observations = Lists.newArrayList();

        Result<Record> result = database.context().resultQuery("select sampleId, svTmb from purity").fetch();
        for (Record record : result) {
            String sampleId = (String) record.getValue(0);
            int svTmb = (Integer) record.getValue(1);

            Sample sample = findSample(samplesToInclude, sampleId);
            if (sample != null) {
                observations.add(ImmutableObservation.builder().type(PercentileType.SV_TMB).sample(sample).value(svTmb).build());
            }
        }

        return observations;
    }

    @Nullable
    private static Sample findSample(@NotNull List<Sample> samples, @NotNull String sampleId) {
        for (Sample sample : samples) {
            if (sample.sampleId().equals(sampleId)) {
                return sample;
            }
        }

        return null;
    }
}

