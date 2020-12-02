package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.CGI;
import static com.hartwig.hmftools.vicc.database.Tables.CGICDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIGDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINDIVIDUALMUTATION;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINFO;
import static com.hartwig.hmftools.vicc.database.Tables.CGIREGION;
import static com.hartwig.hmftools.vicc.database.Tables.CGISTRAND;
import static com.hartwig.hmftools.vicc.database.Tables.CGITRANSCRIPT;

import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class CgiDAOFunctions {

    private CgiDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Cgi cgi) {
        int id = context.insertInto(CGI,
                CGI.GENE,
                CGI.BIOMARKER,
                CGI.ALTERATION,
                CGI.ALTERATIONTYPE,
                CGI.ASSOCIATION,
                CGI.DRUG,
                CGI.DRUGFAMILY,
                CGI.DRUGFULLNAME,
                CGI.DRUGSTATUS,
                CGI.TARGETING,
                CGI.PRIMARYTUMORTYPE,
                CGI.METASTATICTUMORTYPE,
                CGI.EVIDENCELEVEL,
                CGI.SOURCE,
                CGI.CURATOR,
                CGI.ASSAYTYPE,
                CGI.VICCENTRYID)
                .values(cgi.gene(),
                        cgi.biomarker(),
                        cgi.alteration(),
                        cgi.alterationType(),
                        cgi.association(),
                        cgi.drug(),
                        cgi.drugFamily(),
                        cgi.drugFullName(),
                        cgi.drugStatus(),
                        cgi.targeting(),
                        cgi.primaryTumorType(),
                        cgi.metastaticTumorType(),
                        cgi.evidenceLevel(),
                        cgi.source(),
                        cgi.curator(),
                        cgi.assayType(),
                        viccEntryId)
                .returning(CGI.ID)
                .fetchOne()
                .getValue(CGI.ID);

        for (String transcript : cgi.transcripts()) {
            context.insertInto(CGITRANSCRIPT, CGITRANSCRIPT.TRANSCRIPT, CGITRANSCRIPT.CGIID).values(transcript, id).execute();
        }

        for (String individualMutation : cgi.individualMutations()) {
            context.insertInto(CGIINDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.INDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.CGIID)
                    .values(individualMutation, id)
                    .execute();
        }

        for (String gDNA : cgi.gDNA()) {
            context.insertInto(CGIGDNA, CGIGDNA.GDNA, CGIGDNA.CGIID).values(gDNA, id).execute();
        }

        for (String cDNA : cgi.cDNA()) {
            context.insertInto(CGICDNA, CGICDNA.CDNA, CGICDNA.CGIID).values(cDNA, id).execute();
        }

        for (String info : cgi.info()) {
            context.insertInto(CGIINFO, CGIINFO.INFO, CGIINFO.CGIID).values(info, id).execute();
        }

        for (String region : cgi.regions()) {
            context.insertInto(CGIREGION, CGIREGION.REGION, CGIREGION.CGIID).values(region, id).execute();
        }

        for (String strand : cgi.strands()) {
            context.insertInto(CGISTRAND, CGISTRAND.STRAND, CGISTRAND.CGIID).values(strand, id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(CGITRANSCRIPT).execute();
        context.deleteFrom(CGIINDIVIDUALMUTATION).execute();
        context.deleteFrom(CGIGDNA).execute();
        context.deleteFrom(CGICDNA).execute();
        context.deleteFrom(CGIINFO).execute();
        context.deleteFrom(CGIREGION).execute();
        context.deleteFrom(CGISTRAND).execute();

        context.deleteFrom(CGI).execute();
    }
}
