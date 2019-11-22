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
                CGI.TARGETING,
                CGI.SOURCE,
                CGI.PRIMARYTUMORTYPE,
                CGI.DRUGSFULLNAME,
                CGI.CURATOR,
                CGI.DRUGFAMILY,
                CGI.ALTERATION,
                CGI.DRUG,
                CGI.BIOMARKER,
                CGI.DRUGSTATUS,
                CGI.GENE,
                CGI.ASSAYTYPE,
                CGI.ALTERATIONTYPE,
                CGI.EVIDENCELEVEL,
                CGI.ASSOCIATION,
                CGI.METASTATICTUMORTYPE,
                CGI.VICCENTRYID)
                .values(cgi.targeting(),
                        cgi.source(),
                        cgi.primary_tumor_type(),
                        cgi.drugsFullName(),
                        cgi.curator(),
                        cgi.drug_family(),
                        cgi.alteration(),
                        cgi.drug(),
                        cgi.biomarker(),
                        cgi.drug_status(),
                        cgi.gene(),
                        cgi.assay_type(),
                        cgi.alteration_type(),
                        cgi.evidence_level(),
                        cgi.association(),
                        cgi.metastatic_Tumor_Type(),
                        viccEntryId)
                .returning(CGI.ID)
                .fetchOne()
                .getValue(CGI.ID);

        for (String cDNA : cgi.cDNA()) {
            context.insertInto(CGICDNA, CGICDNA.CDNA, CGICDNA.CGIID).values(cDNA, id).execute();
        }

        for (String individualMutation : cgi.individual_mutation()) {
            context.insertInto(CGIINDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.INDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.CGIID)
                    .values(individualMutation, id)
                    .execute();
        }

        for (String gDNA : cgi.gDNA()) {
            context.insertInto(CGIGDNA, CGIGDNA.GDNA, CGIGDNA.CGIID).values(gDNA, id).execute();
        }

        for (String transcript : cgi.transcript()) {
            context.insertInto(CGITRANSCRIPT, CGITRANSCRIPT.TRANSCRIPT, CGITRANSCRIPT.CGIID).values(transcript, id).execute();
        }

        for (String strand : cgi.strand()) {
            context.insertInto(CGISTRAND, CGISTRAND.STRAND, CGISTRAND.CGIID).values(strand, id).execute();
        }

        for (String info : cgi.info()) {
            context.insertInto(CGIINFO, CGIINFO.INFO, CGIINFO.CGIID).values(info, id).execute();
        }

        for (String region : cgi.region()) {
            context.insertInto(CGIREGION, CGIREGION.REGION, CGIREGION.CGIID).values(region, id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(CGICDNA).execute();
        context.deleteFrom(CGIINDIVIDUALMUTATION).execute();
        context.deleteFrom(CGIGDNA).execute();
        context.deleteFrom(CGITRANSCRIPT).execute();
        context.deleteFrom(CGISTRAND).execute();
        context.deleteFrom(CGIINFO).execute();
        context.deleteFrom(CGIREGION).execute();

        context.deleteFrom(CGI).execute();
    }
}
