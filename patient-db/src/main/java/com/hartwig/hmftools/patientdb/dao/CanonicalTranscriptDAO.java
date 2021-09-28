package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CANONICALTRANSCRIPT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;

class CanonicalTranscriptDAO {

    @NotNull
    private final DSLContext context;

    CanonicalTranscriptDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String assembly, @NotNull List<CanonicalTranscript> transcripts) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(CANONICALTRANSCRIPT).where(CANONICALTRANSCRIPT.ASSEMBLY.eq(assembly)).execute();

        for (List<CanonicalTranscript> split : Iterables.partition(transcripts, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep19 inserter = context.insertInto(CANONICALTRANSCRIPT,
                    CANONICALTRANSCRIPT.ASSEMBLY,
                    CANONICALTRANSCRIPT.GENE,
                    CANONICALTRANSCRIPT.GENEID,
                    CANONICALTRANSCRIPT.CHROMOSOMEBAND,
                    CANONICALTRANSCRIPT.CHROMOSOME,
                    CANONICALTRANSCRIPT.GENESTART,
                    CANONICALTRANSCRIPT.GENEEND,
                    CANONICALTRANSCRIPT.TRANSCRIPTID,
                    CANONICALTRANSCRIPT.TRANSCRIPTSTART,
                    CANONICALTRANSCRIPT.TRANSCRIPTEND,
                    CANONICALTRANSCRIPT.EXONS,
                    CANONICALTRANSCRIPT.EXONSTART,
                    CANONICALTRANSCRIPT.EXONEND,
                    CANONICALTRANSCRIPT.EXONBASES,
                    CANONICALTRANSCRIPT.CODINGSTART,
                    CANONICALTRANSCRIPT.CODINGEND,
                    CANONICALTRANSCRIPT.CODINGBASES,
                    CANONICALTRANSCRIPT.STRAND,
                    CANONICALTRANSCRIPT.MODIFIED);
            split.forEach(x -> addRecord(timestamp, inserter, assembly, x));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep19 inserter, @NotNull String assembly,
            @NotNull CanonicalTranscript transcript) {
        inserter.values(assembly,
                transcript.geneName(),
                transcript.geneId(),
                transcript.chromosomeBand(),
                transcript.chromosome(),
                transcript.geneStart(),
                transcript.geneEnd(),
                transcript.transName(),
                transcript.start(),
                transcript.end(),
                transcript.exonCount(),
                transcript.exonStart(),
                transcript.exonEnd(),
                transcript.exonBases(),
                transcript.codingStart(),
                transcript.codingEnd(),
                transcript.codingBases(),
                transcript.strand(),
                timestamp);
    }
}
