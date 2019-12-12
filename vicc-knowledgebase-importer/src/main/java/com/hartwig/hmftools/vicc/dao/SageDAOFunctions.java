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
                SAGE.GENE,
                SAGE.ENTREZID,
                SAGE.CLINICALMANIFESTATION,
                SAGE.RESPONSETYPE,
                SAGE.EVIDENCELABEL,
                SAGE.DRUGLABELS,
                SAGE.GERMLINEORSOMATIC,
                SAGE.PUBLICATIONURL,
                SAGE.VICCENTRYID)
                .values(sage.gene(),
                        sage.entrezId(),
                        sage.clinicalManifestation(),
                        sage.responseType(),
                        sage.evidenceLabel(),
                        sage.drugLabels(),
                        sage.germlineOrSomatic(),
                        sage.publicationUrl(),
                        viccEntryId)
                .execute();
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(SAGE).execute();
    }
}
