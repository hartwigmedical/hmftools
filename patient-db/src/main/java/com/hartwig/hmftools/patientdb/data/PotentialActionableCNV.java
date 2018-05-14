package com.hartwig.hmftools.patientdb.data;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Record;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PotentialActionableCNV {
    @NotNull
    public abstract String sampleId();

    @Nullable
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract CopyNumberAlteration alteration();

    @NotNull
    public static PotentialActionableCNV of(@NotNull final Record mysqlRecord) {
        final int copyNumberValue = (int) Math.max(0, Math.round(mysqlRecord.get(GENECOPYNUMBER.MINCOPYNUMBER)));
        return ImmutablePotentialActionableCNV.of(mysqlRecord.get(GENECOPYNUMBER.SAMPLEID),
                mysqlRecord.get(BASELINE.PRIMARYTUMORLOCATION),
                mysqlRecord.get(GENECOPYNUMBER.GENE),
                CopyNumberAlteration.fromCopyNumber(copyNumberValue));
    }
}
