package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.JAX;
import static com.hartwig.hmftools.vicc.database.Tables.JAXINDICATION;
import static com.hartwig.hmftools.vicc.database.Tables.JAXMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXREFERENCE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTHERAPY;

import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxReference;

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
                JAX.IDJAXENTRY,
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

        context.insertInto(JAXINDICATION, JAXINDICATION.SOURCE, JAXINDICATION.IDINDICATION, JAXINDICATION.NAME, JAXINDICATION.JAXID)
                .values(jax.indication().source(), jax.indication().id(), jax.indication().name(), id)
                .execute();

        for (JaxReference references : jax.references()) {
            context.insertInto(JAXREFERENCE,
                    JAXREFERENCE.URL,
                    JAXREFERENCE.IDREFERENCE,
                    JAXREFERENCE.PUBMEDID,
                    JAXREFERENCE.TITLE,
                    JAXREFERENCE.JAXID).values(references.url(), references.id(), references.pubMedId(), references.title(), id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(JAXMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTHERAPY).execute();
        context.deleteFrom(JAXINDICATION).execute();
        context.deleteFrom(JAXREFERENCE).execute();

        context.deleteFrom(JAX).execute();
    }
}
