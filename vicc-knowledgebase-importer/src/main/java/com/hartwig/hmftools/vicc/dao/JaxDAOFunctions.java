package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.JAX;
import static com.hartwig.hmftools.vicc.database.Tables.JAXINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXREFERENCES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTHERAPY;

import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxReferences;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class JaxDAOFunctions {

    private JaxDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Jax jax) {
        int id = context.insertInto(JAX,
                JAX.RESPONSETYPE,
                JAX.APPROVALSTATUS,
                JAX.EVIDENCETYPE,
                JAX.EFFICACYEVIDENCE,
                JAX.IDJAXSOURCE,
                JAX.VICCENTRYID)
                .values(jax.responseType(), jax.approvalStatus(), jax.evidenceType(), jax.efficacyEvidence(), jax.id(), viccEntryId)
                .returning(JAX.ID)
                .fetchOne()
                .getValue(JAX.ID);

        context.insertInto(JAXMOLECULARPROFILE,
                JAXMOLECULARPROFILE.PROFILENAME,
                JAXMOLECULARPROFILE.IDMOLECULARPROFILE,
                JAXMOLECULARPROFILE.JAXID).values(jax.molecularProfile().profileName(), jax.molecularProfile().id(), id).execute();

        context.insertInto(JAXTHERAPY, JAXTHERAPY.THERAPYNAME, JAXTHERAPY.IDTHERAPY, JAXTHERAPY.JAXID)
                .values(jax.therapy().therapyName(), jax.therapy().id(), id)
                .execute();

        context.insertInto(JAXINDICATIONS, JAXINDICATIONS.SOURCE, JAXINDICATIONS.IDINDICATIONS, JAXINDICATIONS.NAME, JAXINDICATIONS.JAXID)
                .values(jax.indications().source(), jax.indications().id(), jax.indications().name(), id)
                .execute();

        for (JaxReferences references : jax.references()) {
            context.insertInto(JAXREFERENCES,
                    JAXREFERENCES.URL,
                    JAXREFERENCES.IDREFERENCES,
                    JAXREFERENCES.PUBMEDID,
                    JAXREFERENCES.TITLE,
                    JAXREFERENCES.JAXID).values(references.url(), references.id(), references.pubMedId(), references.title(), id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(JAX).execute();
        context.deleteFrom(JAXMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTHERAPY).execute();
        context.deleteFrom(JAXINDICATIONS).execute();
        context.deleteFrom(JAXREFERENCES).execute();
    }
}
