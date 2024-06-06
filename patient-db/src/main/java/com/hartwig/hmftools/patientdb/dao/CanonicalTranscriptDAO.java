package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CANONICALTRANSCRIPT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptUtils;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;

class CanonicalTranscriptDAO {

    private final DSLContext context;

    CanonicalTranscriptDAO(final DSLContext context) {
        this.context = context;
    }

    void write(final String refGenomeVersion, final List<GeneData> geneDataList, final List<TranscriptData> transcripts)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(CANONICALTRANSCRIPT).where(CANONICALTRANSCRIPT.ASSEMBLY.eq(refGenomeVersion)).execute();

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

        for(int i = 0; i < geneDataList.size(); ++i)
        {
            GeneData geneData = geneDataList.get(i);
            TranscriptData transData = transcripts.get(i);

            addRecord(timestamp, inserter, refGenomeVersion, geneData ,transData);
        }

        inserter.execute();
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep19 inserter, final String refGenomeVersion,
            final GeneData geneData, final TranscriptData tranData)
    {
        int exonBases = tranData.exons().stream().mapToInt(x -> x.baseLength()).sum();
        int codingBases = TranscriptUtils.codingBaseLength(tranData);

        inserter.values(refGenomeVersion,
                geneData.GeneName,
                geneData.GeneId,
                geneData.KaryotypeBand,
                geneData.Chromosome,
                geneData.GeneStart,
                geneData.GeneEnd,
                tranData.TransName,
                tranData.TransStart,
                tranData.TransEnd,
                tranData.exons().size(),
                tranData.TransStart,
                tranData.TransEnd,
                exonBases,
                tranData.CodingStart != null ? tranData.CodingStart : 0,
                tranData.CodingEnd != null ? tranData.CodingEnd : 0,
                codingBases,
                geneData.Strand,
                timestamp);
    }
}
