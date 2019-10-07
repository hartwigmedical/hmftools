package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.SAGE;

import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class SageDAOFunctions {

    private SageDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Sage sage) {
        context.insertInto(SAGE,
                SAGE.ENTREZID,
                SAGE.CLINICALMANIFESTATION,
                SAGE.PUBLICATIONURL,
                SAGE.GERMLINEORSOMATIC,
                SAGE.EVIDENCELABEL,
                SAGE.DRUGLABEL,
                SAGE.RESPONSETYPE,
                SAGE.GENE,
                SAGE.VICCENTRYID)
                .values(sage.entrezId(),
                        sage.clinicalManifestation(),
                        sage.publicationUrl(),
                        sage.germlineOrSomatic(),
                        sage.evidenceLabel(),
                        sage.drugLabels(),
                        sage.responseType(),
                        sage.gene(),
                        viccEntryId)
                .execute();
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(SAGE).execute();
    }
}
