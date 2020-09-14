package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.AMBERANONYMOUS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.AMBERMAPPING;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.AMBERPATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.AMBERSAMPLE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.AmberAnonymous;
import com.hartwig.hmftools.common.amber.AmberMapping;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.common.amber.ImmutableAmberAnonymous;
import com.hartwig.hmftools.common.amber.ImmutableAmberMapping;
import com.hartwig.hmftools.common.amber.ImmutableAmberPatient;
import com.hartwig.hmftools.common.amber.ImmutableAmberSample;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep3;
import org.jooq.InsertValuesStep6;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Result;

class AmberDAO {

    @NotNull
    private final DSLContext context;

    AmberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void truncatePatients() {
        context.delete(AMBERPATIENT).execute();
    }

    void truncateMappings() {
        context.delete(AMBERMAPPING).execute();
    }

    @NotNull
    List<AmberPatient> readPatients() {
        final List<AmberPatient> result = Lists.newArrayList();
        final Result<Record> queryResult = context.select().from(AMBERPATIENT).fetch();
        for (Record record : queryResult) {
            result.add(ImmutableAmberPatient.builder()
                    .patientId(record.get(AMBERPATIENT.PATIENTID))
                    .sample(record.get(AMBERPATIENT.SAMPLEID))
                    .build());
        }

        return result;
    }

    @NotNull
    List<AmberSample> readSamples() {
        final List<AmberSample> result = Lists.newArrayList();
        final Result<Record> queryResult = context.select().from(AMBERSAMPLE).fetch();

        for (Record record : queryResult) {
            int size = record.valuesRow().size();
            final byte[] entries = new byte[size - 2];
            for (int i = 0; i < size - 2; i++) {
                entries[i] = record.get(i + 2, Byte.class);
            }

            result.add(ImmutableAmberSample.builder().sampleId(record.get(AMBERSAMPLE.SAMPLEID)).entries(entries).build());
        }

        return result;
    }

    void writePatients(List<AmberPatient> patients) {
        if (patients.isEmpty()) {
            return;
        }

        final Set<String> samples = patients.stream().map(AmberPatient::sample).collect(Collectors.toSet());

        context.delete(AMBERPATIENT).where(AMBERPATIENT.SAMPLEID.in(samples)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());
        InsertValuesStep3 inserter = context.insertInto(AMBERPATIENT, AMBERPATIENT.MODIFIED, AMBERPATIENT.PATIENTID, AMBERPATIENT.SAMPLEID);

        for (AmberPatient amberPatient : patients) {
            inserter.values(timestamp, amberPatient.patientId(), amberPatient.sample());
        }

        inserter.execute();

    }

    void writeAnonymous(List<AmberAnonymous> patients) {
        if (patients.isEmpty()) {
            return;
        }
        context.truncate(AMBERANONYMOUS).execute();
        Timestamp timestamp = new Timestamp(new Date().getTime());
        InsertValuesStep3 inserter = context.insertInto(AMBERANONYMOUS, AMBERANONYMOUS.MODIFIED, AMBERANONYMOUS.SAMPLEID, AMBERANONYMOUS.HMFSAMPLEID);
        for (AmberAnonymous amberPatient : patients) {
            inserter.values(timestamp, amberPatient.sampleId(), amberPatient.hmfSampleId());
        }
        inserter.execute();
    }

    List<AmberAnonymous> readAnonymous() {

        final List<AmberAnonymous> result = Lists.newArrayList();
        final Result<Record> queryResult = context.select().from(AMBERANONYMOUS).fetch();

        for (Record record : queryResult) {
            result.add(ImmutableAmberAnonymous.builder()
                    .sampleId(record.get(AMBERANONYMOUS.SAMPLEID))
                    .hmfSampleId(record.get(AMBERANONYMOUS.HMFSAMPLEID))
                    .build());
        }

        return result;
    }

    List<AmberMapping> readMapping() {
        final List<AmberMapping> result = Lists.newArrayList();
        final Result<Record> queryResult = context.select().from(AMBERMAPPING).fetch();

        for (Record record : queryResult) {
            result.add(ImmutableAmberMapping.builder()
                    .firstSample(record.get(AMBERMAPPING.FIRSTSAMPLEID))
                    .secondSample(record.get(AMBERMAPPING.SECONDSAMPLEID))
                    .matches(record.get(AMBERMAPPING.MATCHES))
                    .sites(record.get(AMBERMAPPING.SITES))
                    .build());
        }

        return result;
    }

    void writeMapping(String sample, List<AmberMapping> mapping) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(AMBERMAPPING).where(AMBERMAPPING.FIRSTSAMPLEID.eq(sample)).execute();
        context.delete(AMBERMAPPING).where(AMBERMAPPING.SECONDSAMPLEID.eq(sample)).execute();

        InsertValuesStep6 inserter = context.insertInto(AMBERMAPPING,
                AMBERMAPPING.MODIFIED,
                AMBERMAPPING.FIRSTSAMPLEID,
                AMBERMAPPING.SECONDSAMPLEID,
                AMBERMAPPING.MATCHES,
                AMBERMAPPING.SITES,
                AMBERMAPPING.LIKELIHOOD);

        for (AmberMapping amberPatient : mapping) {
            inserter.values(timestamp,
                    amberPatient.firstSample(),
                    amberPatient.secondSample(),
                    amberPatient.matches(),
                    amberPatient.sites(),
                    amberPatient.likelihood());
        }

        inserter.execute();
    }

    void writeIdentity(AmberSample identity) {
        byte[] entries = identity.entries();
        if (entries.length != AMBERSAMPLE.fields().length - 2) {
            throw new IllegalArgumentException(
                    "Identity has " + entries.length + " sites but " + (AMBERSAMPLE.fields().length - 2) + " are required");

        }

        context.delete(AMBERSAMPLE).where(AMBERSAMPLE.SAMPLEID.eq(identity.sampleId())).execute();
        InsertValuesStepN inserter = context.insertInto(AMBERSAMPLE,
                AMBERSAMPLE.MODIFIED,
                AMBERSAMPLE.SAMPLEID,
                AMBERSAMPLE.SITE1,
                AMBERSAMPLE.SITE2,
                AMBERSAMPLE.SITE3,
                AMBERSAMPLE.SITE4,
                AMBERSAMPLE.SITE5,
                AMBERSAMPLE.SITE6,
                AMBERSAMPLE.SITE7,
                AMBERSAMPLE.SITE8,
                AMBERSAMPLE.SITE9,
                AMBERSAMPLE.SITE10,
                AMBERSAMPLE.SITE11,
                AMBERSAMPLE.SITE12,
                AMBERSAMPLE.SITE13,
                AMBERSAMPLE.SITE14,
                AMBERSAMPLE.SITE15,
                AMBERSAMPLE.SITE16,
                AMBERSAMPLE.SITE17,
                AMBERSAMPLE.SITE18,
                AMBERSAMPLE.SITE19,
                AMBERSAMPLE.SITE20,
                AMBERSAMPLE.SITE21,
                AMBERSAMPLE.SITE22,
                AMBERSAMPLE.SITE23,
                AMBERSAMPLE.SITE24,
                AMBERSAMPLE.SITE25,
                AMBERSAMPLE.SITE26,
                AMBERSAMPLE.SITE27,
                AMBERSAMPLE.SITE28,
                AMBERSAMPLE.SITE29,
                AMBERSAMPLE.SITE30,
                AMBERSAMPLE.SITE31,
                AMBERSAMPLE.SITE32,
                AMBERSAMPLE.SITE33,
                AMBERSAMPLE.SITE34,
                AMBERSAMPLE.SITE35,
                AMBERSAMPLE.SITE36,
                AMBERSAMPLE.SITE37,
                AMBERSAMPLE.SITE38,
                AMBERSAMPLE.SITE39,
                AMBERSAMPLE.SITE40,
                AMBERSAMPLE.SITE41,
                AMBERSAMPLE.SITE42,
                AMBERSAMPLE.SITE43,
                AMBERSAMPLE.SITE44,
                AMBERSAMPLE.SITE45,
                AMBERSAMPLE.SITE46,
                AMBERSAMPLE.SITE47,
                AMBERSAMPLE.SITE48,
                AMBERSAMPLE.SITE49,
                AMBERSAMPLE.SITE50,
                AMBERSAMPLE.SITE51,
                AMBERSAMPLE.SITE52,
                AMBERSAMPLE.SITE53,
                AMBERSAMPLE.SITE54,
                AMBERSAMPLE.SITE55,
                AMBERSAMPLE.SITE56,
                AMBERSAMPLE.SITE57,
                AMBERSAMPLE.SITE58,
                AMBERSAMPLE.SITE59,
                AMBERSAMPLE.SITE60,
                AMBERSAMPLE.SITE61,
                AMBERSAMPLE.SITE62,
                AMBERSAMPLE.SITE63,
                AMBERSAMPLE.SITE64,
                AMBERSAMPLE.SITE65,
                AMBERSAMPLE.SITE66,
                AMBERSAMPLE.SITE67,
                AMBERSAMPLE.SITE68,
                AMBERSAMPLE.SITE69,
                AMBERSAMPLE.SITE70,
                AMBERSAMPLE.SITE71,
                AMBERSAMPLE.SITE72,
                AMBERSAMPLE.SITE73,
                AMBERSAMPLE.SITE74,
                AMBERSAMPLE.SITE75,
                AMBERSAMPLE.SITE76,
                AMBERSAMPLE.SITE77,
                AMBERSAMPLE.SITE78,
                AMBERSAMPLE.SITE79,
                AMBERSAMPLE.SITE80,
                AMBERSAMPLE.SITE81,
                AMBERSAMPLE.SITE82,
                AMBERSAMPLE.SITE83,
                AMBERSAMPLE.SITE84,
                AMBERSAMPLE.SITE85,
                AMBERSAMPLE.SITE86,
                AMBERSAMPLE.SITE87,
                AMBERSAMPLE.SITE88,
                AMBERSAMPLE.SITE89,
                AMBERSAMPLE.SITE90,
                AMBERSAMPLE.SITE91,
                AMBERSAMPLE.SITE92,
                AMBERSAMPLE.SITE93,
                AMBERSAMPLE.SITE94,
                AMBERSAMPLE.SITE95,
                AMBERSAMPLE.SITE96,
                AMBERSAMPLE.SITE97,
                AMBERSAMPLE.SITE98,
                AMBERSAMPLE.SITE99,
                AMBERSAMPLE.SITE100);

        List<Object> collection = Lists.newArrayList();
        collection.add(new Timestamp(new Date().getTime()));
        collection.add(identity.sampleId());
        for (final byte entry : entries) {
            collection.add(entry);
        }

        inserter.values(collection);
        inserter.execute();
    }


    void deleteAmberRecordsForSample(@NotNull String sample) {
        context.delete(AMBERPATIENT).where(AMBERPATIENT.SAMPLEID.eq(sample)).execute();
        context.delete(AMBERMAPPING).where(AMBERMAPPING.FIRSTSAMPLEID.eq(sample)).execute();
        context.delete(AMBERMAPPING).where(AMBERMAPPING.SECONDSAMPLEID.eq(sample)).execute();
        context.delete(AMBERSAMPLE).where(AMBERSAMPLE.SAMPLEID.eq(sample)).execute();
    }
}