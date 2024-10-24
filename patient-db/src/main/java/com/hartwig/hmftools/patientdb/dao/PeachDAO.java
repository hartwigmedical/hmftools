package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PEACHGENOTYPE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;
import org.jooq.Record;
import org.jooq.Result;

class PeachDAO
{
    private final DSLContext context;

    PeachDAO(final DSLContext context)
    {
        this.context = context;
    }

    public List<PeachGenotype> readPeachGenotypes(final String sample)
    {
        Result<Record> result = context.select().from(PEACHGENOTYPE).where(PEACHGENOTYPE.SAMPLEID.eq(sample)).fetch();
        List<PeachGenotype> genotypes = Lists.newArrayList();
        for(Record record : result)
        {
            PeachGenotype genotype = ImmutablePeachGenotype.builder()
                    .gene(record.getValue(PEACHGENOTYPE.GENE))
                    .allele(record.getValue(PEACHGENOTYPE.HAPLOTYPE))
                    .alleleCount(record.getValue(PEACHGENOTYPE.COUNT))
                    .function(record.getValue(PEACHGENOTYPE.FUNCTION))
                    .linkedDrugs(record.getValue(PEACHGENOTYPE.LINKEDDRUGS))
                    .urlPrescriptionInfo(record.getValue(PEACHGENOTYPE.URLPRESCRIPTIONINFO))
                    .build();
            genotypes.add(genotype);
        }
        return genotypes;
    }

    void writePeach(final String sample, final List<PeachGenotype> peachGenotypes)
    {
        deletePeachForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for(List<PeachGenotype> genotypes : Iterables.partition(peachGenotypes, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep8 inserter = context.insertInto(PEACHGENOTYPE,
                    PEACHGENOTYPE.MODIFIED,
                    PEACHGENOTYPE.SAMPLEID,
                    PEACHGENOTYPE.GENE,
                    PEACHGENOTYPE.HAPLOTYPE,
                    PEACHGENOTYPE.COUNT,
                    PEACHGENOTYPE.FUNCTION,
                    PEACHGENOTYPE.LINKEDDRUGS,
                    PEACHGENOTYPE.URLPRESCRIPTIONINFO);
            genotypes.forEach(x -> addGenotype(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addGenotype(
            final Timestamp timestamp, final InsertValuesStep8 inserter, final String sample, final PeachGenotype genotype)
    {
        inserter.values(timestamp,
                sample,
                genotype.gene(),
                genotype.allele(),
                genotype.alleleCount(),
                genotype.function(),
                genotype.linkedDrugs(),
                genotype.urlPrescriptionInfo());
    }

    void deletePeachForSample(@NotNull String sample)
    {
        context.delete(PEACHGENOTYPE).where(PEACHGENOTYPE.SAMPLEID.eq(sample)).execute();
    }
}
